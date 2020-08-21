#ifndef STATICCHECKPOINTING
#define STATICCHECKPOINTING

#include "petsc.h"
#include "topOpt/WangCheckpointing.hh"
#include <limits>

class StaticCheckpointing : public WangCheckpointing {

public:
  StaticCheckpointing(PetscInt,PetscInt,Vec);
  ~StaticCheckpointing();

  PetscErrorCode makeCheckpoint(PetscInt,Vec,void*) override;
  PetscErrorCode reset(void*) override;
private:
  struct CheckpointInfo {
    PetscInt step, level;
  };
  PetscInt numTimestep;
  PetscInt curCheck;
  CheckpointInfo* forwardCheckpoints;
};

#endif
