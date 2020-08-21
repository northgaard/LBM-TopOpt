#ifndef THERMALLBSOLVER2D
#define THERMALLBSOLVER2D

#include "petsc.h"
#include "NewLBSolver2d.hh"
#include "core/definitions2d.hh"
#include "core/macros.hh"
#include "core/meta.hh"
#include "core/codiheader.hh"
#include "dynamics/thermalLBLoops2d.hh"
#include "adjointDynamics/adjointThermalLBLoops2d.hh"

#include <type_traits>

template <class ThermalCollisionOperator>
class ThermalLBSolver2d : public NewLBSolver2d {

public:
  PetscErrorCode make(MPI_Comm _comm, const NewLBSolverInfo2d& info, ThermalCollisionOperator _op)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    communicator = _comm;
    colOp = _op;
    ierr = createPetscObjects(info,IsoLattice::numDOF+ThermalLattice::numDOF,4,
                              ThermalCollisionOperator::numAdditionalFields);
    CHKERRQ(ierr);
    ierr = createBoundingBoxes(info); CHKERRQ(ierr);

    ierr = DMDASetFieldName(macroGrid,0,"rho"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,1,"ux"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,2,"uy"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,3,"T"); CHKERRQ(ierr);

    getCollisionLoops();
    streamBySwappingLoop =
      new (std::nothrow) ThermalStreamBySwappingLoop2d<ThermalCollisionOperator>
      (localBoundingBoxGhosted); CHKNEWPTR(streamBySwappingLoop);
    initializeAtEquilibriumLoop =
      new (std::nothrow) ThermalEquilibriumInitializationLoop2d<ThermalCollisionOperator>
      (localBoundingBox); CHKNEWPTR(initializeAtEquilibriumLoop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode uniformMacroInitialization(const ThermalMacros2d&);
  PetscErrorCode getInitialMacroArray(const Box2d&,Box2d*,ThermalMacros2d***);
  PetscErrorCode restoreInitialMacroArray(const Box2d&,Box2d*,ThermalMacros2d***);

  template <class T = ThermalCollisionOperator>
  typename std::enable_if<T::numAdditionalFields,PetscErrorCode>::type
  makeCodiAdjoint(AdjointLBSolver* adjSolver)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = commonAdjointSetup(adjSolver); CHKERRQ(ierr);

    auto adjColOp = colOp.getCodiAdjoint();
    using AdjointCollisionLoop = CodiAdjointObstThermalCollisionLoop2d<T,decltype(adjColOp)>;
    AdjointCollisionLoop* adjColLoop =
      new (std::nothrow) AdjointCollisionLoop(localBoundingBox,adjColOp);
    CHKNEWPTR(adjColLoop);
    ierr = setAdjointCollisionLoop(adjColLoop,adjSolver); CHKERRQ(ierr);

    using AdjointStreamingLoop = NewAdjointThermalStreamingLoop2d<T>;
    AdjointStreamingLoop* adjStreamLoop =
      new (std::nothrow) AdjointStreamingLoop(localBoundingBox,localBoundingBoxGhosted);
    CHKNEWPTR(adjStreamLoop);
    ierr = setAdjointStreamingLoop(adjStreamLoop,adjSolver); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  template <class T = ThermalCollisionOperator>
  typename std::enable_if<T::numAdditionalFields,PetscErrorCode>::type
  makeSourceAdjoint(AdjointLBSolver* adjSolver)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = commonAdjointSetup(adjSolver); CHKERRQ(ierr);

    auto adjColOp = colOp.getSourceAdjoint();
    using AdjointCollisionLoop = SourceAdjointObstThermalCollisionLoop2d<T,decltype(adjColOp)>;
    AdjointCollisionLoop* adjColLoop =
      new (std::nothrow) AdjointCollisionLoop(localBoundingBox,adjColOp);
    CHKNEWPTR(adjColLoop);
    ierr = setAdjointCollisionLoop(adjColLoop,adjSolver); CHKERRQ(ierr);

    using AdjointStreamingLoop = NewAdjointThermalStreamingLoop2d<T>;
    AdjointStreamingLoop* adjStreamLoop =
      new (std::nothrow) AdjointStreamingLoop(localBoundingBox,localBoundingBoxGhosted);
    CHKNEWPTR(adjStreamLoop);
    ierr = setAdjointStreamingLoop(adjStreamLoop,adjSolver); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ThermalLBSolver2d(){}
private:
  template <class T = ThermalCollisionOperator>
  typename std::enable_if<T::numAdditionalFields,PetscErrorCode>::type
  getCollisionLoops()
  {
    PetscFunctionBeginUser;
    collideAndSwapLoop =
      new (std::nothrow) ObstThermalCollisionAndSwapLoop2d<ThermalCollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndSwapLoop);
    collideAndStreamLoop =
      new (std::nothrow) ObstThermalCollideAndStreamSingleLoop2d<ThermalCollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndStreamLoop);
    PetscFunctionReturn(0);
  }

  // This should be extracted to common functionality for LB solvers
  PetscErrorCode commonAdjointSetup(AdjointLBSolver* adjSolver)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    // Ensure design grid has been created
    ierr = getDesignGrid(nullptr); CHKERRQ(ierr);
    ierr = setUpAdjointSolver(adjSolver); CHKERRQ(ierr);
    for (const auto& boundary : boundaryContainer){
      AdjointLBBoundaryLoop* adjBoundaryLoop;
      ierr = boundary->getAdjointBoundaryLoop(&adjBoundaryLoop); CHKERRQ(ierr);
      ierr = addAdjointBoundaryCondition(adjBoundaryLoop,adjSolver); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  ThermalCollisionOperator colOp;
};

template <class ThermalCollisionOperator>
PetscErrorCode ThermalLBSolver2d<ThermalCollisionOperator>::
uniformMacroInitialization(const ThermalMacros2d& uinit)
{
  PetscErrorCode ierr;
  ThermalMacros2d** macInit;
  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,&macInit); CHKERRQ(ierr);

  for (auto jj : localBoundingBox.yRange){
    for (auto ii : localBoundingBox.xRange){
      macInit[jj][ii].rho = uinit.rho;
      macInit[jj][ii].ux = uinit.ux;
      macInit[jj][ii].uy = uinit.uy;
      macInit[jj][ii].T = uinit.T;
    }
  }

  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,&macInit); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

template <class ThermalCollisionOperator>
PetscErrorCode ThermalLBSolver2d<ThermalCollisionOperator>::
getInitialMacroArray(const Box2d& inputBox,
                     Box2d* outLocalBox,
                     ThermalMacros2d*** inpPtr)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  boxIntersection(inputBox,localBoundingBox,*outLocalBox);
  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,inpPtr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

template <class ThermalCollisionOperator>
PetscErrorCode ThermalLBSolver2d<ThermalCollisionOperator>::
restoreInitialMacroArray(const Box2d&,Box2d*,
                         ThermalMacros2d*** inpPtr)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,inpPtr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif
