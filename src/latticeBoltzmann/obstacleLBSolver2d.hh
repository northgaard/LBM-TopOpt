#ifndef OBSTACLELBSOLVER2D
#define OBSTACLELBSOLVER2D

#include "petsc.h"
#include "core/definitions2d.hh"
#include "core/geometry2d.hh"
#include "descriptors/collisionDescriptorBases2d.hh"
#include "latticeBoltzmann/obstacleLBSolverBase2d.hh"
#include <memory>

class ObstacleLBSolverInfo2d {

  friend class ObstacleLBSolver2d;

public:

  ObstacleLBSolverInfo2d(){
    obstacleDOF = 1;
    boundaryX = DM_BOUNDARY_NONE;
    boundaryY = DM_BOUNDARY_NONE;
  }

  PetscInt nx,ny;
  PetscInt latticeDOF;
  PetscInt macroDOF;

private:

  PetscInt obstacleDOF;
  PetscDMBoundary boundaryX, boundaryY;

};

class ObstacleLBSolver2d : public ObstacleLBSolverBase2d {

public:

  ObstacleLBSolver2d() : ObstacleLBSolverBase2d(){}
  PetscErrorCode initializeSolver(const ObstacleLBSolverInfo2d&);

  inline void setCollisionFunction(ObstacleDynamicsFunction _fun)
  {
    collisionFunc = _fun;
  }

  PetscErrorCode streamAndCollide() override;
  PetscErrorCode collideAndStream() override;
  PetscErrorCode collide() override;

private:

  ObstacleDynamicsFunction collisionFunc;

};

#endif
