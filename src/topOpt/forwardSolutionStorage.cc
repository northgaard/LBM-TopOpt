#include "forwardSolutionStorage.hh"

PetscErrorCode ForwardSolutionStorage::getForwardSolution(PetscInt index,
							  Vec& vec)
{

  vec = forwardSolutions[index];
  return 0;

}

PetscErrorCode ForwardSolutionStorage::saveForwardSolution(PetscInt index,
							   Vec vec)
{

  PetscErrorCode ierr;
  ierr = VecCopy(vec,forwardSolutions[index]); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ForwardSolutionStorage::allocate(PetscInt _nV, Vec repVec)
{

  PetscErrorCode ierr;
  numVecs = _nV;
  ierr = VecDuplicateVecs(repVec,numVecs,&forwardSolutions); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ForwardSolutionStorage::deallocate()
{

  PetscErrorCode ierr;
  ierr = VecDestroyVecs(numVecs,&forwardSolutions); CHKERRQ(ierr);

  return 0;

}

