#ifndef LBLOOP
#define LBLOOP

#include "petsc.h"

class LBLoop {

public:
  virtual void execute(void*,void*,void*) = 0;
  virtual ~LBLoop(){}
};

#endif
