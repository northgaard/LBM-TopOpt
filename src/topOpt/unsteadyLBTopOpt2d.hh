#ifndef UNSTEADYLBTOPOPT2D
#define UNSTEADYLBTOPOPT2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "core/memoryHandling.hh"
#include "latticeBoltzmann/obstacleLBSolver2d.hh"
#include "adjointLatticeBoltzmann/adjointLBSolver2d.hh"
#include "objectiveFunctions/objectiveFunction2d.hh"
#include "topOpt/LBTopOpt2d.hh"
#include "topOpt/forwardSolutionStorage.hh"
#include "topOpt/filter.hh"
#include <memory>

class UnsteadyLBTopOpt2d : public LBTopOpt2d {

  using Objective = std::unique_ptr<ObjectiveFunction2d>;

public:

  UnsteadyLBTopOpt2d(ObstacleLBSolver2d& _solv, AdjointLBSolver2d& _asolv,
		     const Objective&  _obj, PetscInt _nt,
		     std::unique_ptr<Filter> _fil = nullptr)
    : LBTopOpt2d(_solv,_asolv,std::move(_fil)), obj(_obj), store(0), numTimesteps(_nt)
  {
    adjointSolver.setNumTimesteps(_nt);
  }

  PetscErrorCode computeObjective();
  PetscErrorCode forwardAndAdjoint();
  PetscErrorCode allocateMemoryForForwardSolution(PetscInt numCheckpoints = -1);
  PetscErrorCode deallocateMemoryForForwardSolution();
  
  PetscScalar getObjective(){ return objectiveValue; }
  

private:

  PetscErrorCode computeObjectiveWithCheckpoints();

  const Objective& obj;
  ForwardSolutionStorage* store;
  PetscInt numTimesteps;
  
};

#endif
