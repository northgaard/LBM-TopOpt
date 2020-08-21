#ifndef ISOTHERMALLBSOLVER2D
#define ISOTHERMALLBSOLVER2D

#include "petsc.h"
#include "NewLBSolver2d.hh"
#include "core/definitions2d.hh"
#include "core/macros.hh"
#include "core/meta.hh"
#include "core/codiheader.hh"
#include "dynamics/isothermalLBLoops2d.hh"
#include "adjointDynamics/adjointIsothermalLBLoops2d.hh"
#include <type_traits>

template <class CollisionOperator>
class IsothermalLBSolver2d : public NewLBSolver2d {

public:

  /* Setup of main solver */
  PetscErrorCode make(MPI_Comm _comm, const NewLBSolverInfo2d& info, CollisionOperator _op)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    communicator = _comm;
    colOp = _op;
    ierr = createPetscObjects(info,Lattice::numDOF,3,CollisionOperator::numAdditionalFields);
    CHKERRQ(ierr);
    ierr = createBoundingBoxes(info); CHKERRQ(ierr);

    ierr = DMDASetFieldName(macroGrid,0,"rho"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,1,"ux"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,2,"uy"); CHKERRQ(ierr);

    ierr = getCollisionLoops(); CHKERRQ(ierr);
    // streamBySwappingLoop =
    //   new (std::nothrow) IsoStreamBySwappingLoop2d<Lattice>(localBoundingBoxGhosted);
    // CHKNEWPTR(streamBySwappingLoop);
    streamLoop =
      new (std::nothrow) IsoStreamLoop2d<Lattice>(localBoundingBox,globalBoundingBox);
    initializeAtEquilibriumLoop =
      new (std::nothrow) IsoEquilibriumInitializationLoop2d
      <typename CollisionOperator::Equilibrium>(localBoundingBox);
    CHKNEWPTR(initializeAtEquilibriumLoop);
    PetscFunctionReturn(0);
  }
  /* Setup of CoDiPack adjoint solver */
  template <class T = CollisionOperator>
  typename std::enable_if<(has_codi_adjoint<T>::value && T::numAdditionalFields),
                          PetscErrorCode>::type
  makeCodiAdjoint(AdjointLBSolver* adjSolver)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = commonAdjointSetup(adjSolver); CHKERRQ(ierr);

    auto adjColOp = colOp.getCodiAdjoint();
    using AdjointCollisionLoop = CodiAdjointObstIsoCollisionLoop2d<T,decltype(adjColOp)>;
    AdjointCollisionLoop* adjColLoop =
      new (std::nothrow) AdjointCollisionLoop(localBoundingBox,adjColOp);
    CHKNEWPTR(adjColLoop);
    ierr = setAdjointCollisionLoop(adjColLoop,adjSolver); CHKERRQ(ierr);

    using AdjointStreamingLoop = AdjointIsoStreamingLoop2d<T>;
    AdjointStreamingLoop* adjStreamLoop =
      new (std::nothrow) AdjointStreamingLoop(localBoundingBox,globalBoundingBox);
    CHKNEWPTR(adjStreamLoop);
    ierr = setAdjointStreamingLoop(adjStreamLoop,adjSolver); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template <class T = CollisionOperator>
  typename std::enable_if<!(has_codi_adjoint<T>::value && T::numAdditionalFields),
                          PetscErrorCode>::type
  makeCodiAdjoint(AdjointLBSolver*)
  {
    static_assert(has_codi_adjoint<T>::value,"The collision operator does not define the getCodiAdjoint() method.\n");
    static_assert(T::numAdditionalFields,"Adjoint LBM is only available for collision operators with an additional field (i.e. porosity for partial bounceback).\n");
    return 0;
  }
  template <class T = CollisionOperator>
  typename std::enable_if<(has_source_adjoint<T>::value && T::numAdditionalFields),
                          PetscErrorCode>::type
  makeSourceAdjoint(AdjointLBSolver* adjSolver)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = commonAdjointSetup(adjSolver); CHKERRQ(ierr);

    auto adjColOp = colOp.getSourceAdjoint();
    using AdjointCollisionLoop = SourceAdjointObstIsoCollisionLoop2d<T,decltype(adjColOp)>;
    AdjointCollisionLoop* adjColLoop =
      new (std::nothrow) AdjointCollisionLoop(localBoundingBox,adjColOp);
    CHKNEWPTR(adjColLoop);
    ierr = setAdjointCollisionLoop(adjColLoop,adjSolver); CHKERRQ(ierr);

    using AdjointStreamingLoop = AdjointIsoStreamingLoop2d<T>;
    AdjointStreamingLoop* adjStreamLoop =
      new (std::nothrow) AdjointStreamingLoop(localBoundingBox,globalBoundingBox);
    CHKNEWPTR(adjStreamLoop);
    ierr = setAdjointStreamingLoop(adjStreamLoop,adjSolver); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template <class T = CollisionOperator>
  typename std::enable_if<!(has_source_adjoint<T>::value && T::numAdditionalFields),
                          PetscErrorCode>::type
  makeSourceAdjoint(AdjointLBSolver* adjSolver)
  {
    static_assert(has_source_adjoint<T>::value,"The collision operator does not define the getSourceAdjoint() method.\n");
    static_assert(T::numAdditionalFields,"Adjoint LBM is only available for collision operators with an additional field (i.e. porosity for partial bounceback).\n");
  }

  /* Initialization of obstacle field */
  // template <class T = CollisionOperator>
  // typename std::enable_if<T::numAdditionalFields,PetscErrorCode>::type
  // getFieldArray(const Box2d& inputBox, Box2d* outLocalBox, PetscScalar**** inpPtr)
  // {
  //   PetscErrorCode ierr;
  //   PetscFunctionBeginUser;
  //   boxIntersection(inputBox,localBoundingBox,*outLocalBox);
  //   ierr = DMDAVecGetArrayDOF(materialGrid,materialLocal,inpPtr); CHKERRQ(ierr);
  //   PetscFunctionReturn(0);
  // }
  // template <class T = CollisionOperator>
  // typename std::enable_if<T::numAdditionalFields,PetscErrorCode>::type
  // restoreFieldArray(const Box2d&,Box2d*,PetscScalar**** inpPtr)
  // {
  //   PetscErrorCode ierr;
  //   PetscFunctionBeginUser;
  //   ierr = DMDAVecRestoreArrayDOF(materialGrid,materialLocal,inpPtr); CHKERRQ(ierr);
  //   ierr = DMLocalToLocalBegin(materialGrid,materialLocal,INSERT_VALUES,
  //                              materialLocal); CHKERRQ(ierr);
  //   ierr = DMLocalToLocalBegin(materialGrid,materialLocal,INSERT_VALUES,
  //                              materialLocal); CHKERRQ(ierr);
  //   PetscFunctionReturn(0);
  // }

  PetscErrorCode uniformMacroInitialization(const IsothermalMacros2d&);
  PetscErrorCode getInitialMacroArray(const Box2d&,Box2d*,IsothermalMacros2d***);
  PetscErrorCode restoreInitialMacroArray(const Box2d&,Box2d*,IsothermalMacros2d***);

  IsothermalLBSolver2d(){}
private:
  using Lattice = typename CollisionOperator::LatticeType;
  CollisionOperator colOp;

  PetscErrorCode commonAdjointSetup(AdjointLBSolver*);
  /* Method for creating collision loops */
  template <class T = CollisionOperator>
  typename std::enable_if<T::numAdditionalFields,PetscErrorCode>::type
  getCollisionLoops()
  {
    PetscFunctionBeginUser;
    collideLoop = new (std::nothrow) ObstIsoCollisionLoop2d<CollisionOperator>
      (localBoundingBox,colOp); CHKNEWPTR(collideLoop);
    // collideAndSwapLoop =
    //   new (std::nothrow) ObstIsoCollisionAndSwapLoop2d<Lattice,CollisionOperator>
    //   (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndSwapLoop);
    // collideAndStreamLoop =
    //   new (std::nothrow) ObstIsoCollideAndStreamSingleLoop2d<Lattice,CollisionOperator>
    //   (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndStreamLoop);
    PetscFunctionReturn(0);
  }
  template <class T = CollisionOperator>
  typename std::enable_if<!T::numAdditionalFields,PetscErrorCode>::type
  getCollisionLoops()
  {
    PetscFunctionBeginUser;
    collideLoop = new (std::nothrow) IsoCollisionLoop2d<CollisionOperator>
      (localBoundingBox,colOp); CHKNEWPTR(collideLoop);
    // collideAndSwapLoop =
    //   new (std::nothrow) IsoCollisionAndSwapLoop2d<Lattice,CollisionOperator>
    //   (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndSwapLoop);
    // collideAndStreamLoop =
    //   new (std::nothrow) IsoCollideAndStreamSingleLoop2d<Lattice,CollisionOperator>
    //   (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndStreamLoop);
    PetscFunctionReturn(0);
  }
};

template <class CollisionOperator>
PetscErrorCode IsothermalLBSolver2d<CollisionOperator>::
uniformMacroInitialization(const IsothermalMacros2d& uinit)
{
  PetscErrorCode ierr;
  IsothermalMacros2d** macInit;
  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,&macInit); CHKERRQ(ierr);

  for (auto jj : localBoundingBox.yRange){
    for (auto ii : localBoundingBox.xRange){
      macInit[jj][ii].rho = uinit.rho;
      macInit[jj][ii].ux = uinit.ux;
      macInit[jj][ii].uy = uinit.uy;
    }
  }

  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,&macInit); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

template <class CollisionOperator>
PetscErrorCode IsothermalLBSolver2d<CollisionOperator>::
getInitialMacroArray(const Box2d& inputBox,
                     Box2d* outLocalBox,
                     IsothermalMacros2d*** inpPtr)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  boxIntersection(inputBox,localBoundingBox,*outLocalBox);
  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,inpPtr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

template <class CollisionOperator>
PetscErrorCode IsothermalLBSolver2d<CollisionOperator>::
restoreInitialMacroArray(const Box2d&,Box2d*,
                         IsothermalMacros2d*** inpPtr)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,inpPtr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

template <class CollisionOperator>
PetscErrorCode IsothermalLBSolver2d<CollisionOperator>::
commonAdjointSetup(AdjointLBSolver* adjSolver)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  // Ensure design grid has been created
  ierr = getDesignGrid(nullptr); CHKERRQ(ierr);
  ierr = setUpAdjointSolver(adjSolver); CHKERRQ(ierr);
  for (const auto& boundary : boundaryContainer){
    AdjointLBBoundaryLoop* adjBoundaryLoop;
    ierr = boundary->getAdjointBoundaryLoop(&adjBoundaryLoop); CHKERRQ(ierr);
    if (adjBoundaryLoop){
      ierr = addAdjointBoundaryCondition(adjBoundaryLoop,adjSolver); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#endif
