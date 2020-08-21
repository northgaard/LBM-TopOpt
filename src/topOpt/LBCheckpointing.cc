#include "LBCheckpointing.hh"

// Static
LBStaticCheckpointing::LBStaticCheckpointing(PetscInt numCheckpoints,
                                             PetscInt numTimestep,
                                             LBSolver& _solv)
  : StaticCheckpointing(numCheckpoints,numTimestep,_solv.distributionsGlobal),
    solver(_solv)
{}

PetscErrorCode LBStaticCheckpointing::advance(void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = solver.collideAndStream(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBStaticCheckpointing::setSolverState(PetscInt curStep, Vec newDist, void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = solver.setDistributions(newDist); CHKERRQ(ierr);
  ierr = solver.setCurrentTimestep(curStep); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBStaticCheckpointing::copyFromSolverState(Vec copy, void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecCopy(solver.distributionsGlobal,copy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Dynamic
LBDynamicCheckpointing::LBDynamicCheckpointing(PetscInt numCheckpoints,
                                               LBSolver& _solv)
  : DynamicCheckpointing(numCheckpoints,_solv.distributionsGlobal), solver(_solv)
{}

PetscErrorCode LBDynamicCheckpointing::advance(void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = solver.collideAndStream(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBDynamicCheckpointing::setSolverState(PetscInt curStep, Vec newDist, void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = solver.setDistributions(newDist); CHKERRQ(ierr);
  ierr = solver.setCurrentTimestep(curStep); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBDynamicCheckpointing::copyFromSolverState(Vec copy, void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecCopy(solver.distributionsGlobal,copy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
