#ifndef WANGCHECKPOINTING
#define WANGCHECKPOINTING

#include "petsc.h"
#include "topOpt/checkpointingBase.hh"
#include <list>
#include <stack>
#include <limits>

class WangCheckpointing : public CheckpointingBase {

protected:
  struct Checkpoint {
    Checkpoint() : step(0), level(0), data(nullptr) {}
    Checkpoint(PetscInt _s, PetscInt _l, Vec _v)
      : step(_s), level(_l), data(_v) {}

    PetscInt step,level;
    Vec data;
  };
public:
  WangCheckpointing(PetscInt,Vec);
  virtual ~WangCheckpointing();

  PetscErrorCode reset(void*) override;
  PetscErrorCode getForwardSolution(PetscInt,Vec,void*) override;
protected:
  std::list<Checkpoint>::iterator computeNewCheckpoints(PetscInt);

  std::list<Checkpoint> checkpoints;
  std::stack<Vec> availableStorage;
  Vec* storage;
  const PetscInt numCheck;
};

#endif
