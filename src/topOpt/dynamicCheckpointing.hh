#ifndef DYNAMICCHECKPOINTING
#define DYNAMICCHECKPOINTING

#include "petsc.h"
#include "topOpt/WangCheckpointing.hh"
#include <limits>

class DynamicCheckpointing : public WangCheckpointing {

public:
  DynamicCheckpointing(PetscInt,Vec);
  ~DynamicCheckpointing();

  PetscErrorCode makeCheckpoint(PetscInt,Vec,void*) override;
private:
  std::list<Checkpoint>::iterator reader;
};

#endif
