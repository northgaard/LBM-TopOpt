#ifndef NAIVECHECKPOINTING
#define NAIVECHECKPOINTING

#include "petsc.h"
#include "topOpt/checkpointingBase.hh"

class NaiveCheckpointing : public CheckpointingBase {

public:

  NaiveCheckpointing(PetscInt,Vec);
  ~NaiveCheckpointing(); 
  PetscErrorCode makeCheckpoint(PetscInt,Vec,void*) override;
  PetscErrorCode getForwardSolution(PetscInt,Vec,void*) override;
  PetscErrorCode reset(void*) override;
private:
  PetscErrorCode advance(void*) override;
  PetscErrorCode setSolverState(PetscInt,Vec,void*) override;
  PetscErrorCode copyFromSolverState(Vec,void*) override;
  Vec* forwardSolutions;
  PetscInt numTimestep;
};

#endif
