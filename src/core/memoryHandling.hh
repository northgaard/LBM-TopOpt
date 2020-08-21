#ifndef MEMORYHANDLING
#define MEMORYHANDLING

#include "petsc.h"

struct _pList {
  void* pointer;
  struct _pList* next;
};

typedef struct _pList* MemoryList;

PetscErrorCode registerPointerForDeallocation(void*);
PetscErrorCode deallocateAllRegisteredPointers();

#endif
