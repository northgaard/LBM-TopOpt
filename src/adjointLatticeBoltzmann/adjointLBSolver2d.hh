#ifndef ADJOINTLBSOLVER2D
#define ADJOINTLBSOLVER2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/definitions2d.hh"
#include "core/geometry2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "objectiveFunctions/objectiveFunction2d.hh"
#include "latticeBoltzmann/obstacleLBSolver2d.hh"
#include <functional>
#include <vector>
#include <memory>

class AdjointLBSolver2d {

protected:

  using BaseDynamicsFunction = std::function<void(PetscScalar***,
						  PetscScalar***)>;
  using ObstacleDynamicsFunction = std::function<void(PetscScalar***,
						      PetscScalar***,
						      PetscScalar**,
						      PetscScalar**)>;
  using AdjointBoundaryFunction = std::function<void(PetscInt,PetscScalar***,
						     PetscScalar***)>;

public:

  AdjointLBSolver2d() : adjointGlobal(0), adjointLocal(0), latticeGrid(0),
			macroGrid(0), obstacleGrid(0), numTimesteps(1) {}

  PetscErrorCode initializeFromForwardSolver(const ObstacleLBSolver2d& solver);

  PetscErrorCode adjointCollide(Vec,Vec,Vec,ObjectiveFunction2d*);
  PetscErrorCode adjointStream();
  
  PetscErrorCode resetAdjoints()
  {

    PetscErrorCode ierr;
    ierr = VecSet(adjointGlobal,0.); CHKERRQ(ierr);
    ierr = VecSet(adjointLocal,0.); CHKERRQ(ierr);
    curTimestep = numTimesteps;

    return 0;

  }

  inline void setAdjointStreamingFunction(BaseDynamicsFunction _fun)
  {
    adjointStreamFunc = _fun;
  }
  inline void setAdjointCollisionFunction(ObstacleDynamicsFunction _fun)
  {
    adjointCollisionFunc = _fun;
  }
  void setNumTimesteps(PetscInt nt){ numTimesteps = nt; }

  PetscErrorCode addBoundaryCondition(const Box2d,
				      const std::unique_ptr<BoundaryDescriptor2d>&);


protected:

  PetscErrorCode createAdjointVectors();
  void boundaryComputations(PetscInt,PetscScalar***,PetscScalar***);

  BaseDynamicsFunction adjointStreamFunc;
  ObstacleDynamicsFunction adjointCollisionFunc;
  std::vector<AdjointBoundaryFunction> adjointBoundaryContainer;

  Vec adjointGlobal;
  Vec adjointLocal;
  DM latticeGrid;
  DM macroGrid;
  DM obstacleGrid;
  Box2d globalBoundingBox;
  Box2d localBoundingBox;
  Box2d localBoundingBoxGhosted;
  PetscInt numTimesteps,curTimestep;

};

PetscErrorCode
addBoundaryCondition(const Box2d,
		     const std::unique_ptr<BoundaryDescriptor2d>&,
		     ObstacleLBSolver2d&,AdjointLBSolver2d&);

#endif
