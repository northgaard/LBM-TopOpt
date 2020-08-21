#ifndef MATERIALFUNCTIONAL
#define MATERIALFUNCTIONAL

#include "petsc.h"

class MaterialFunctional {

public:
  virtual void computeFunctional(Vec,PetscScalar*) const = 0;
  virtual void computeSensitivities(Vec) const = 0;
};

#endif
