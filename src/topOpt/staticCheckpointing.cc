#include "staticCheckpointing.hh"

StaticCheckpointing::StaticCheckpointing(PetscInt _nc, PetscInt _nt, Vec _vec)
  : WangCheckpointing(_nc,_vec), numTimestep(_nt)
{

  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject) storage[0], &comm);
  ierr = PetscMalloc1(numCheck,&forwardCheckpoints);
  CHKERRABORT(comm,ierr);

  computeNewCheckpoints(numTimestep-1);
  PetscInt id = 0;

  for (auto itr = checkpoints.rbegin(); itr != checkpoints.rend(); ++itr){
    forwardCheckpoints[id].step = (*itr).step;
    forwardCheckpoints[id].level = (*itr).level;
    ++id;
  }

  reset(nullptr);

}

StaticCheckpointing::~StaticCheckpointing()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject) storage[0], &comm);
  ierr = PetscFree(forwardCheckpoints);
  CHKERRABORT(comm,ierr);
}

PetscErrorCode StaticCheckpointing::makeCheckpoint(PetscInt index, Vec check, void*)
{
  PetscFunctionBeginUser;
  if (forwardCheckpoints[curCheck].step == index){
    PetscErrorCode ierr;
    ierr = VecCopy(check,availableStorage.top()); CHKERRQ(ierr);
    checkpoints.emplace_front(forwardCheckpoints[curCheck].step,
                              forwardCheckpoints[curCheck].level,
                              availableStorage.top());
    availableStorage.pop();
    ++curCheck;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode StaticCheckpointing::reset(void*)
{
  WangCheckpointing::reset(nullptr);
  curCheck = 0;
  return 0;
}
