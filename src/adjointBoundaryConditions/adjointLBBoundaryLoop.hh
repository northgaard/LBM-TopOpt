#ifndef ADJOINTLBBOUNDARYLOOP
#define ADJOINTLBBOUNDARYLOOP

#include "petsc.h"

class AdjointLBBoundaryLoop {

public:
  virtual void execute(PetscInt,void*,void*) const = 0;
  virtual ~AdjointLBBoundaryLoop(){}
};

#endif
