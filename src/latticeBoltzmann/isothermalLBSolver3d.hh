#ifndef ISOTHERMALLBSOLVER3D
#define ISOTHERMALLBSOLVER3D

#include "petsc.h"
#include "LBSolver3d.hh"
#include "core/definitions3d.hh"
#include "core/macros.hh"
#include "core/codiheader.hh"
#include "dynamics/isothermalLBLoops3d.hh"
#include "traits/collisionTraits3d.hh"
#include <type_traits>

template <class CollisionOperator>
class IsothermalLBSolver3d : public LBSolver3d {

  using Lattice = typename get_lattice<CollisionOperator>::type;
public:
  PetscErrorCode make(MPI_Comm _comm, const LBSolverInfo3d& info, CollisionOperator _op)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    communicator = _comm;
    colOp = _op;
    ierr = createPetscObjects(info,Lattice::numDOF,4,
                              get_num_additional_fields<CollisionOperator>::value);
    CHKERRQ(ierr);
    ierr = createBoundingBoxes(info); CHKERRQ(ierr);

    ierr = DMDASetFieldName(macroGrid,0,"rho"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,1,"ux"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,2,"uy"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(macroGrid,3,"uz"); CHKERRQ(ierr);

    ierr = getCollisionLoops(); CHKERRQ(ierr);
    streamBySwappingLoop =
      new (std::nothrow) IsoStreamBySwappingLoop3d<CollisionOperator>
      (localBoundingBoxGhosted); CHKNEWPTR(streamBySwappingLoop);
    initializeAtEquilibriumLoop =
      new (std::nothrow) IsoEquilibriumInitializationLoop3d<CollisionOperator>
      (localBoundingBox); CHKNEWPTR(initializeAtEquilibriumLoop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode uniformMacroInitialization(const IsothermalMacros3d& uinit)
  {
    PetscErrorCode ierr;
    IsothermalMacros3d*** macInit;
    PetscFunctionBeginUser;
    ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,&macInit); CHKERRQ(ierr);

    for (auto kk : localBoundingBox.zRange){
      for (auto jj : localBoundingBox.yRange){
        for (auto ii : localBoundingBox.xRange){
          macInit[kk][jj][ii].rho = uinit.rho;
          macInit[kk][jj][ii].ux = uinit.ux;
          macInit[kk][jj][ii].uy = uinit.uy;
          macInit[kk][jj][ii].uz = uinit.uz;
        }
      }
    }

    ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,&macInit); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getInitialMacroArray(const Box3d& inputBox,
                                      Box3d* outLocalBox,
                                      IsothermalMacros3d**** inpPtr)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    boxIntersection(inputBox,localBoundingBox,*outLocalBox);
    ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,inpPtr); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode restoreInitialMacroArray(const Box3d&,Box3d*,
                                          IsothermalMacros3d**** inpPtr)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,inpPtr); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  IsothermalLBSolver3d(){}
private:
  template <class T = CollisionOperator>
  typename std::enable_if<get_num_additional_fields<T>::value,PetscErrorCode>::type
  getCollisionLoops()
  {
    PetscFunctionBeginUser;
    collideLoop = new (std::nothrow) ObstIsoCollisionLoop3d<CollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideLoop);
    collideAndSwapLoop =
      new (std::nothrow) ObstIsoCollisionAndSwapLoop3d<CollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndSwapLoop);
    collideAndStreamLoop =
      new (std::nothrow) ObstIsoCollideAndStreamSingleLoop3d<CollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndStreamLoop);
    PetscFunctionReturn(0);
  }
  template <class T = CollisionOperator>
  typename std::enable_if<!get_num_additional_fields<T>::value,PetscErrorCode>::type
  getCollisionLoops()
  {
    PetscFunctionBeginUser;
    collideLoop = new (std::nothrow) IsoCollisionLoop3d<CollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideLoop);
    collideAndSwapLoop =
      new (std::nothrow) IsoCollisionAndSwapLoop3d<CollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndSwapLoop);
    collideAndStreamLoop =
      new (std::nothrow) IsoCollideAndStreamSingleLoop3d<CollisionOperator>
      (localBoundingBoxGhosted,colOp); CHKNEWPTR(collideAndStreamLoop);
    PetscFunctionReturn(0);
  }
  CollisionOperator colOp;
};

#endif
