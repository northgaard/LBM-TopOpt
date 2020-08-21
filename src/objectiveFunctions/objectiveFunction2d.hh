#ifndef OBJECTIVEFUNCTION2D
#define OBJECTIVEFUNCTION2D

#include "petsc.h"

class ObjectiveFunction2d {

public:

  virtual void
  evaluate(PetscScalar***,PetscScalar&) const = 0;
  virtual void
  adjointStateSource(PetscInt,PetscScalar***,PetscScalar***) const = 0;
  virtual void
  adjointDesignSource(PetscInt,PetscScalar**,PetscScalar**) const {}

};

#endif
