#ifndef CHECKPOINTINGBASE
#define CHECKPOINTINGBASE

#include "petsc.h"
#include "latticeBoltzmann/LBSolverBase.hh"

class CheckpointingBase {

public:
  virtual ~CheckpointingBase(){};
  virtual PetscErrorCode makeCheckpoint(PetscInt,Vec,void* userContext) = 0;
  virtual PetscErrorCode getForwardSolution(PetscInt,Vec,void* userContext) = 0;
  virtual PetscErrorCode reset(void* userContext) = 0;
protected:
  virtual PetscErrorCode advance(void* userContext) = 0;
  virtual PetscErrorCode setSolverState(PetscInt,Vec,void* userContext) = 0;
  virtual PetscErrorCode copyFromSolverState(Vec,void*) = 0;
};

#endif
