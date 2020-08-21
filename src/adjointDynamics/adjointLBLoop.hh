#ifndef ADJOINTLBLOOP
#define ADJOINTLBLOOP

#include "petsc.h"

class AdjointLBLoop {

public:
  virtual void execute(void*,void*,void*,void*) = 0;
  virtual void executeNoSensitivities(void*,void*,void*){}
  virtual ~AdjointLBLoop(){}
};

#endif
