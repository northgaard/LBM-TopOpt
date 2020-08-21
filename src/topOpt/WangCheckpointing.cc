#include "WangCheckpointing.hh"

WangCheckpointing::WangCheckpointing(PetscInt _nc, Vec _vec)
  : numCheck(_nc)
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject) _vec, &comm);
  ierr = VecDuplicateVecs(_vec,numCheck,&storage);
  CHKERRABORT(comm,ierr);

  for (PetscInt ii = 0; ii < numCheck; ++ii){
    availableStorage.push(storage[ii]);
  }
}

WangCheckpointing::~WangCheckpointing()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject) storage[0], &comm);
  ierr = VecDestroyVecs(numCheck,&storage);
  CHKERRABORT(comm,ierr);
}

PetscErrorCode WangCheckpointing::reset(void*)
{
  for (const auto& chk : checkpoints){
    availableStorage.push(chk.data);
  }
  checkpoints.clear();

  return 0;
}

PetscErrorCode WangCheckpointing::getForwardSolution(PetscInt index,
                                                     Vec vector, void* userContext)
{
  std::list<Checkpoint>::iterator init;
  init = checkpoints.begin();
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  if ((*init).step != index){
    init = computeNewCheckpoints(index);

    // PetscPrintf(PETSC_COMM_WORLD,"Checkpoint at init: %d %d\n",(*init).step,(*init).level);
    // PetscPrintf(PETSC_COMM_WORLD,"New checkpoints:\n");
    // for (auto itr = checkpoints.rbegin(); itr != checkpoints.rend(); ++itr){
    //   PetscPrintf(PETSC_COMM_WORLD,"%d %d\n",(*itr).step,(*itr).level);
    // }

    PetscInt curStep = (*init).step;
    // PetscPrintf(PETSC_COMM_WORLD,"Starting recomputation from step %d\n",curStep);
    // solver->setDistributions((*init).data);
    ierr = setSolverState(curStep,(*init).data,userContext); CHKERRQ(ierr);
    --init;

    for (;;){
      ierr = advance(userContext); CHKERRQ(ierr);
      ++curStep;
      // PetscPrintf(PETSC_COMM_WORLD,"Advanced to step %d\n",curStep);
      if ((*init).step == curStep){
        // PetscPrintf(PETSC_COMM_WORLD,"Placing new checkpoint at step %d\n",curStep);
        if ((*init).step == index){
          break;
        }
        ierr = copyFromSolverState((*init).data,userContext); CHKERRQ(ierr);
        --init;
      }
    }
    ierr = copyFromSolverState(vector,userContext); CHKERRQ(ierr);
  } else {
    ierr = VecCopy((*init).data,vector);
  }

  availableStorage.push((*init).data);
  checkpoints.erase(init);
  // PetscPrintf(PETSC_COMM_WORLD,"Grapping checkpoint at step %d\n",index);
  PetscFunctionReturn(0);
}

std::list<WangCheckpointing::Checkpoint>::iterator
WangCheckpointing::computeNewCheckpoints(PetscInt index)
{
  static const PetscInt maxLevel =
    std::numeric_limits<PetscInt>::max();

  if (checkpoints.empty()){
    checkpoints.emplace_front(0,maxLevel,availableStorage.top());
    availableStorage.pop();
  }

  std::list<Checkpoint>::iterator startingPoint = checkpoints.begin();
  std::list<Checkpoint>::iterator reader = checkpoints.begin();
  std::list<Checkpoint>::iterator movable;

  PetscInt numSteps = index - (*startingPoint).step;
  PetscInt initStep = (*startingPoint).step;
  PetscInt newLevel;

  for (PetscInt ii = initStep; ii < initStep+numSteps; ++ii){

    if (availableStorage.size()){
      checkpoints.emplace_front(ii+1,0,availableStorage.top());
      availableStorage.pop();
      reader = checkpoints.begin();
      continue;
    }

    movable = reader;
    ++movable;

    if ((*movable).level < (*reader).level){
      checkpoints.emplace_front(ii+1,0,(*movable).data);
      checkpoints.erase(movable);
    } else {
      reader = checkpoints.begin();
      newLevel = (*reader).level + 1;
      checkpoints.emplace_front(ii+1,newLevel,(*reader).data);
      checkpoints.erase(reader);
      reader = checkpoints.begin();
    }
  }
  return startingPoint;
}
