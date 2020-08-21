#ifndef BINARYIO
#define BINARYIO

#include "petsc.h"
#include <string>

PetscErrorCode writeVectorToBinary(const std::string&,const std::string&,Vec,
                                   MPI_Comm = PETSC_COMM_WORLD);
PetscErrorCode loadVectorFromBinary(const std::string&,const std::string&,Vec,
                                    MPI_Comm = PETSC_COMM_WORLD);

#endif
