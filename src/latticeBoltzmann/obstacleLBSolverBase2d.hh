#ifndef OBSTACLELBSOLVERBASE2D
#define OBSTACLELBSOLVERBASE2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"
#include <functional>
#include <string>

class ObstacleLBSolverBase2d : public LBSolverBase2d {

public:

  Vec obstacleGlobal;
  DM obstacleGrid;

protected:

  using ObstacleDynamicsFunction = std::function<void(PetscScalar***,
						      PetscScalar***,
						      PetscScalar**)>;

  ObstacleLBSolverBase2d() : LBSolverBase2d(),
			     obstacleGlobal(0), obstacleGrid(0) {}
  PetscErrorCode createObstacleObjects(PetscInt);
  
public:

  PetscErrorCode addGeometricObstacle(const Box2d&, PetscScalar);
  PetscErrorCode outputDomain(std::string);
  PetscErrorCode setDomain(Vec);

};

#endif
