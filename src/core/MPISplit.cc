#include "MPISplit.hh"
#include <utility>

/*
  Implementation file for MPIEvenSplit.
  Author: Sebastian Arlund Nørgaard (sebnorg@login.gbar.dtu.dk)
*/

/*
  Constructor
  Input:
      _sp -- integer equal to the number of subcommunicators
*/
MPIEvenSplit::MPIEvenSplit(PetscMPIInt _sp) : split(_sp), setUpScatterCalled(false)
{
  PetscMPIInt* numProcsPerComm = new PetscMPIInt[split];
  PetscMPIInt worldSize;
  PetscMPIInt worldRank;
  MPI_Comm_size(PETSC_COMM_WORLD,&worldSize);
  MPI_Comm_rank(PETSC_COMM_WORLD,&worldRank);
  PetscMPIInt div = worldSize / split;
  PetscMPIInt remain = worldSize % split;

  // Determine how to split the processes
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    numProcsPerComm[ii] = div;
  }
  for (PetscMPIInt ii = 0; ii < remain; ++ii){
    ++numProcsPerComm[ii];
  }
  std::pair<PetscMPIInt,PetscMPIInt>* grouping =
    new std::pair<PetscMPIInt,PetscMPIInt>[split];
  PetscMPIInt beg = 0;
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    grouping[ii].first = beg;
    grouping[ii].second = grouping[ii].first + numProcsPerComm[ii] - 1;
    beg += numProcsPerComm[ii];
  }
  // Determine the ID of the group of processors
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    if (grouping[ii].first <= worldRank && worldRank <= grouping[ii].second){
      splitID = ii;
      break;
    }
  }

  // Create communicator
  MPI_Comm_split(PETSC_COMM_WORLD,splitID,worldRank,&activeCommunicator);

  // Calculate world roots
  worldRoots = new PetscMPIInt[split];
  worldRoots[0] = 0;
  for (PetscMPIInt ii = 1; ii < split; ++ii){
    worldRoots[ii] = worldRoots[ii-1] + numProcsPerComm[ii-1];
  }

  delete [] numProcsPerComm;
  delete [] grouping;
}

/* Destructor */
MPIEvenSplit::~MPIEvenSplit()
{
  MPI_Comm_free(&activeCommunicator);
  delete [] worldRoots;
  if (setUpScatterCalled){
    PetscErrorCode ierr;
    for (PetscMPIInt ii = 0; ii < split; ++ii){
      ierr = VecDestroy(&(dummiesSplit[ii]));
      CHKERRABORT(PETSC_COMM_WORLD,ierr);
      ierr = VecScatterDestroy(&(scatters[ii]));
      CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
    delete [] dummiesSplit;
    delete [] scatters;
  }
}

/*
  Sets up PETSc objects to handle scattering between the subset communicator and "the world"
  (e.g. PETSC_COMM_WORLD).
  Input:
      splitVec -- Vector created on the subset communicator
      worldVec -- Vector create on PETSC_COMM_WORLD
  Both input vectors should be "the same", i.e. have the same size etc. They should only
  differ in their communicators.
*/
PetscErrorCode MPIEvenSplit::setUpScatter(Vec splitVec, Vec worldVec)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  if (setUpScatterCalled){
    PetscFunctionReturn(0);
  }
  // Create dummy vectors
  dummiesSplit = new Vec[split];
  scatters = new VecScatter[split];
  PetscInt splitTotalSize, splitLocalSize, splitBlockSize;
  /*
    Are there potential problems with an uneven split of processors,
    since their vectors can then have different block sizes?
  */
  ierr = VecGetSize(splitVec,&splitTotalSize); CHKERRQ(ierr);
  ierr = VecGetLocalSize(splitVec,&splitLocalSize); CHKERRQ(ierr);
  ierr = VecGetBlockSize(splitVec,&splitBlockSize); CHKERRQ(ierr);
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    if (ii == splitID){
      ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,splitBlockSize,splitLocalSize,
                                   splitTotalSize,NULL,&(dummiesSplit[ii]));
      CHKERRQ(ierr);
    } else {
      ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,splitBlockSize,0,
                                   splitTotalSize,NULL,&(dummiesSplit[ii]));
      CHKERRQ(ierr);
    }
  }
  // Create scatters
  scatters = new VecScatter[split];
  PetscInt low,high;
  ierr = VecGetOwnershipRange(splitVec,&low,&high); CHKERRQ(ierr);
  IS* indexSets = new IS[split];
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    if (ii == splitID){
      ierr = ISCreateStride(PETSC_COMM_WORLD,high-low,low,1,&(indexSets[ii]));
      CHKERRQ(ierr);
    } else {
      ierr = ISCreateStride(PETSC_COMM_WORLD,0,0,0,&(indexSets[ii]));
      CHKERRQ(ierr);
    }
    ierr = VecScatterCreate(dummiesSplit[ii],indexSets[ii],worldVec,
                            indexSets[ii],&(scatters[ii])); CHKERRQ(ierr);
  }
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    ierr = ISDestroy(&(indexSets[ii])); CHKERRQ(ierr);
  }
  delete [] indexSets;
  setUpScatterCalled = true;
  PetscFunctionReturn(0);
}

/*
  Scatters from local communicator to an array of world vectors
  Input:
      splitVec -- Vector created on subset communicator
      worldVecs -- Array of vectors on PETSC_COMM_WORLD. Length must be greater or equal
                   to the number of subcommunicators (_sp in constructor)
  The vectors must have the same layout as those given to setUpScatter(), the scatter is
  from splitVec ---> worldVecs[splitID]
*/
PetscErrorCode MPIEvenSplit::scatterToWorld(Vec splitVec, Vec* worldVecs)
{
  PetscErrorCode ierr;
  PetscScalar* array;
  PetscFunctionBeginUser;
  if (!setUpScatterCalled){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"You must call setUpScatter before you can call scatterToWorld.\n");
  }
  // Set data in dummy vector
  ierr = VecGetArray(splitVec,&array); CHKERRQ(ierr);
  ierr = VecPlaceArray(dummiesSplit[splitID],array); CHKERRQ(ierr);
  // Perform scatters to world vectors
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    ierr = VecScatterBegin(scatters[ii],dummiesSplit[ii],worldVecs[ii],INSERT_VALUES,
                           SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatters[ii],dummiesSplit[ii],worldVecs[ii],INSERT_VALUES,
                         SCATTER_FORWARD); CHKERRQ(ierr);
  }
  ierr = VecResetArray(dummiesSplit[splitID]); CHKERRQ(ierr);
  ierr = VecRestoreArray(splitVec,&array); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Similar to scatterToWorld() above, but only scatters from one subset communicator to
  a world vector. Essentially
  if (IDToScatter == splitID){
      splitVec ---> worldVec
  }
  This is useful if for example you have a volume constraint on a specific realization,
  the requirements on the Vecs are otherwise identical to scatterToWorld() above
*/
PetscErrorCode MPIEvenSplit::scatterSingleToWorld(Vec splitVec, Vec worldVec,
                                                  PetscMPIInt IDToScatter)
{
  PetscErrorCode ierr;
  PetscScalar* array;
  PetscFunctionBeginUser;
  if (!setUpScatterCalled){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"You must call setUpScatter before you can call scatterSingleToWorld.\n");
  }
  if (IDToScatter == splitID){
    ierr = VecGetArray(splitVec,&array); CHKERRQ(ierr);
    ierr = VecPlaceArray(dummiesSplit[IDToScatter],array); CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(scatters[IDToScatter],dummiesSplit[IDToScatter],worldVec,
                         INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(scatters[IDToScatter],dummiesSplit[IDToScatter],worldVec,
                       INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  if (IDToScatter == splitID){
    ierr = VecResetArray(dummiesSplit[IDToScatter]); CHKERRQ(ierr);
    ierr = VecRestoreArray(splitVec,&array); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  Scatter a vector from the world to all subset communicators. This allow you to scatter your
  updated design from MMA to all realizations
  Input:
      worldVec -- The vector on PETSC_COMM_WORLD
      splitVec -- The subset communicators copy of the world vector
*/
PetscErrorCode MPIEvenSplit::scatterCommonFromWorld(Vec worldVec, Vec splitVec)
{
  PetscErrorCode ierr;
  PetscScalar* array;
  PetscFunctionBeginUser;
  if (!setUpScatterCalled){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"You must call setUpScatter before you can call scatterCommonFromWorld.\n");
  }
  // Set data in dummy vector
  ierr = VecGetArray(splitVec,&array); CHKERRQ(ierr);
  ierr = VecPlaceArray(dummiesSplit[splitID],array); CHKERRQ(ierr);
  // Perform scatters from world vector
  for (PetscMPIInt ii = 0; ii < split; ++ii){
    ierr = VecScatterBegin(scatters[ii],worldVec,dummiesSplit[ii],INSERT_VALUES,
                           SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatters[ii],worldVec,dummiesSplit[ii],INSERT_VALUES,
                         SCATTER_REVERSE); CHKERRQ(ierr);
  }
  ierr = VecResetArray(dummiesSplit[splitID]); CHKERRQ(ierr);
  ierr = VecRestoreArray(splitVec,&array); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Gets necessary information from the object.
  Output:
      theComm -- The subset communicator
           id -- Integer identifying the subset communicator, essentially a rank for the group
                 of processors.
        roots -- Array of length equal to the number of subset communicators, roots[x] contains
                 the world rank of the process which has rank 0 in the subset communicator with
                 id == x
*/
void MPIEvenSplit::getCommunicationInfo(MPI_Comm* theComm, PetscMPIInt* id, PetscMPIInt** roots)
{
  if (theComm){
    *theComm = activeCommunicator;
  }
  if (id){
    *id = splitID;
  }
  if (roots){
    *roots = worldRoots;
  }
}
