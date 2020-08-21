#ifndef COLLISIONDESCRIPTORBASES2D
#define COLLISIONDESCRIPTORBASES2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include <functional>

class CollisionDescriptorBase2d {

public:

  using BaseDynamicsFunction = std::function<void(PetscScalar***,PetscScalar***)>;

  virtual PetscErrorCode
  getInitializationLoop(Box2d,BaseDynamicsFunction&) const = 0;
  virtual PetscErrorCode
  getStreamingLoop(Box2d,Box2d,BaseDynamicsFunction&) const = 0;
  virtual PetscErrorCode
  getCollisionLoop(Box2d,BaseDynamicsFunction&) const = 0;
  virtual PetscInt getDistributionDOF() const = 0 ;
  virtual PetscInt getMacroDOF() const = 0;

};

class ObstacleCollisionDescriptorBase2d {

public:

  using BaseDynamicsFunction = std::function<void(PetscScalar***,PetscScalar***)>;
  using ObstacleDynamicsFunction =
    std::function<void(PetscScalar***,PetscScalar***,PetscScalar**)>;
  using AdjointObstacleDynamicsFunction =
    std::function<void(PetscScalar***,PetscScalar***,
		       PetscScalar**,PetscScalar**)>;

  virtual PetscErrorCode
  getInitializationLoop(Box2d,BaseDynamicsFunction&) const = 0;
  virtual PetscErrorCode
  getStreamingLoop(Box2d,Box2d,BaseDynamicsFunction&) const = 0;
  virtual PetscErrorCode
  getAdjointStreamingLoop(Box2d,Box2d,BaseDynamicsFunction&) const;
  virtual PetscErrorCode
  getCollisionLoop(Box2d,ObstacleDynamicsFunction&) const = 0;
  virtual PetscErrorCode
  getAdjointCollisionLoop(Box2d,AdjointObstacleDynamicsFunction&) const;
  virtual PetscInt getDistributionDOF() const = 0;
  virtual PetscInt getMacroDOF() const = 0;

};

#endif
