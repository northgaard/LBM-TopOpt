#include "petsc.h"

#define CHKNEWPTR(p) do {if (PetscUnlikely(!(p))) return PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME_CXX,__FILE__,PETSC_ERR_MEM,PETSC_ERROR_INITIAL,"Operator new was unable to allocate requested memory.\n");} while (0)
