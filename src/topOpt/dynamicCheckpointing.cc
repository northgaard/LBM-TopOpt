#include "dynamicCheckpointing.hh"

DynamicCheckpointing::DynamicCheckpointing(PetscInt _nc, Vec _vec)
  : WangCheckpointing(_nc,_vec)
{}

DynamicCheckpointing::~DynamicCheckpointing(){}

PetscErrorCode DynamicCheckpointing::makeCheckpoint(PetscInt index, Vec check,
                                                    void* userContext)
{
  static const PetscInt maxLevel =
    std::numeric_limits<PetscInt>::max();
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  // Check for initial checkpoint
  if (checkpoints.empty()){
    // TODO: Add error if index != 0
    checkpoints.emplace_front(index,maxLevel,availableStorage.top());
    reader = checkpoints.begin();
    // ierr = copyFromSolverState(availableStorage.top(),userContext); CHKERRQ(ierr);
    ierr = VecCopy(check,availableStorage.top()); CHKERRQ(ierr);
    availableStorage.pop();
    // PetscPrintf(PETSC_COMM_WORLD,"Placing 0 checkpoint.\n");
    PetscFunctionReturn(0);
  }
  // Place checkpoint if not yet filled up
  if (availableStorage.size()){
    checkpoints.emplace_front(index,0,availableStorage.top());
    reader = checkpoints.begin();
    ierr = VecCopy(check,availableStorage.top()); CHKERRQ(ierr);
    availableStorage.pop();
    PetscFunctionReturn(0);
  }
  // Else we have to redistribute
  std::list<Checkpoint>::iterator movable = reader;
  ++movable;

  if ((*movable).level < (*reader).level){
    checkpoints.emplace_front(index,0,(*movable).data);
    ierr = VecCopy(check,(*movable).data); CHKERRQ(ierr);
    checkpoints.erase(movable);
  } else {
    reader = checkpoints.begin();
    PetscInt newLevel = (*reader).level + 1;
    checkpoints.emplace_front(index,newLevel,(*reader).data);
    ierr = VecCopy(check,(*reader).data); CHKERRQ(ierr);
    checkpoints.erase(reader);
    reader = checkpoints.begin();
  }

  PetscFunctionReturn(0);
}
