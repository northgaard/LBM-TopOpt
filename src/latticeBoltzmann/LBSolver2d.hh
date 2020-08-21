#ifndef LBSOLVER2D
#define LBSOLVER2D

#include "petsc.h"
#include "core/definitions2d.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"
#include <memory>

class LBSolverInfo2d {

  friend class LBSolver2d;

public:

  LBSolverInfo2d(){
    boundaryX = DM_BOUNDARY_NONE;
    boundaryY = DM_BOUNDARY_NONE;
  }
  
  PetscInt nx,ny;
  PetscInt latticeDOF;
  PetscInt macroDOF;

private:

  PetscDMBoundary boundaryX, boundaryY;

};

class LBSolver2d : public LBSolverBase2d {

public:

  LBSolver2d() : LBSolverBase2d(){}
  ~LBSolver2d(){}
  PetscErrorCode initializeSolver(const LBSolverInfo2d&);

  inline void setCollisionFunction(BaseDynamicsFunction _fun)
  {
    collisionFunc = _fun;
  }

  PetscErrorCode streamAndCollide() override;
  PetscErrorCode collideAndStream() override;
  PetscErrorCode collide() override;

private:

  BaseDynamicsFunction collisionFunc;
    
};

#endif
