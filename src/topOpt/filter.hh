#ifndef FILTER
#define FILTER

#include "petsc.h"
#include <memory>

class Filter {

public:

  virtual PetscErrorCode filterDesign(Vec,Vec) = 0;
  virtual PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) = 0;
  virtual PetscErrorCode outputIntermediateFields(PetscInt)
  {
    return 0;
  }

  virtual ~Filter(){}

};

#endif
