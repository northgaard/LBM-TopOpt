#include "petsc.h"

/*
  Class for splitting computational resources in evenly sized communicators,
  useful for e.g. robust optimization where many realizations of the same problem
  can be run in an embarrassingly parallel fashion.

  See implementation file for documentation.

  Author: Sebastian Arlund Nørgaard (sebnorg@mek.dtu.dk)
*/
class MPIEvenSplit {

public:
  MPIEvenSplit(PetscMPIInt _sp);
  ~MPIEvenSplit();

  PetscErrorCode setUpScatter(Vec,Vec);
  PetscErrorCode scatterToWorld(Vec,Vec*);
  PetscErrorCode scatterSingleToWorld(Vec,Vec,PetscMPIInt);
  PetscErrorCode scatterCommonFromWorld(Vec,Vec);
  void getCommunicationInfo(MPI_Comm*,PetscMPIInt*,PetscMPIInt**);
private:
  MPI_Comm activeCommunicator;
  Vec* dummiesSplit;
  VecScatter* scatters;
  PetscMPIInt splitID;
  PetscMPIInt split;
  PetscMPIInt* worldRoots;
  bool setUpScatterCalled;
};
