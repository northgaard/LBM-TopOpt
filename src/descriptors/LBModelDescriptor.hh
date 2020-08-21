#ifndef LBMODELDESCRIPTOR
#define LBMODELDESCRIPTOR

#include "petsc.h"
#include <functional>

class LBModelDescriptor {

  using DynamicsFunction = std::function<void(void*,void*,void*)>;

public:

  PetscInt virtual getNumDOF() = 0;
  PetscInt virtual getNumMacros() = 0;
  PetscInt virtual getNumAdditionalFields() = 0;

  PetscErrorCode virtual getCollideAndSwapFunc(DynamicsFunction*) = 0;
  PetscErrorCode virtual getCollideAndStreamFunc(DynamicsFunction*) = 0;
  // PetscErrorCode virtual getCollideFunc(DynamicsFunction*) = 0;
  PetscErrorCode virtual getStreamBySwappingFunc(DynamicsFunction*) = 0;
  // PetscErrorCode virtual getStreamFunc(DynamicsFunction*) = 0;
  // PetscErrorCode virtual getComputeMacrosFunc(DynamicsFunction*) = 0;
  PetscErrorCode virtual getInitializeAtEquilibriumFunc(DynamicsFunction*) = 0;
  
};

#endif
