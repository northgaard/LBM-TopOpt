#include "naiveCheckpointing.hh"

NaiveCheckpointing::NaiveCheckpointing(PetscInt _nt, Vec repVec) : numTimestep(_nt)
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject) repVec, &comm);
  ierr = VecDuplicateVecs(repVec,numTimestep,&forwardSolutions);
  CHKERRABORT(comm,ierr);
}

NaiveCheckpointing::~NaiveCheckpointing()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject) forwardSolutions[0], &comm);
  ierr = VecDestroyVecs(numTimestep,&forwardSolutions);
  CHKERRABORT(comm,ierr);
}

PetscErrorCode NaiveCheckpointing::makeCheckpoint(PetscInt index,Vec check,void*)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecCopy(check,forwardSolutions[index]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NaiveCheckpointing::getForwardSolution(PetscInt index,
                                                      Vec check,void*)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  ierr = VecCopy(forwardSolutions[index],check); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NaiveCheckpointing::reset(void*)
{
  return 0;
}

// Naive checkpointing doesn't need a solver
PetscErrorCode NaiveCheckpointing::advance(void*)
{
  return 0;
}

PetscErrorCode NaiveCheckpointing::setSolverState(PetscInt,Vec,void*)
{
  return 0;
}

PetscErrorCode NaiveCheckpointing::copyFromSolverState(Vec,void*)
{
  return 0;
}
