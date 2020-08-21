#ifndef FORWARDSTORAGE
#define FORWARDSTORAGE

#include "petsc.h"

class ForwardSolutionStorage {

public:

  PetscErrorCode getForwardSolution(PetscInt,Vec&);
  PetscErrorCode saveForwardSolution(PetscInt,Vec);
  PetscErrorCode allocate(PetscInt,Vec);
  PetscErrorCode deallocate();

private:

  Vec* forwardSolutions;
  PetscInt numVecs;

};

#endif
