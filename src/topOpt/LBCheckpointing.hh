#ifndef LBSTATICCHECKPOINTING
#define LBSTATICCHECKPOINTING

#include "petsc.h"
#include "topOpt/staticCheckpointing.hh"
#include "topOpt/dynamicCheckpointing.hh"
#include "latticeBoltzmann/LBSolver.hh"

class LBStaticCheckpointing : public StaticCheckpointing {

public:
  LBStaticCheckpointing(PetscInt,PetscInt,LBSolver&);
private:
  PetscErrorCode advance(void*) override;
  PetscErrorCode setSolverState(PetscInt,Vec,void*) override;
  PetscErrorCode copyFromSolverState(Vec,void*) override;
  LBSolver& solver;
};

class LBDynamicCheckpointing : public DynamicCheckpointing {

public:
  LBDynamicCheckpointing(PetscInt,LBSolver&);
private:
  PetscErrorCode advance(void*) override;
  PetscErrorCode setSolverState(PetscInt,Vec,void*) override;
  PetscErrorCode copyFromSolverState(Vec,void*) override;
  LBSolver& solver;
};

#endif
