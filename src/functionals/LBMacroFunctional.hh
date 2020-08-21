#ifndef LBMACROFUNCTIONAL
#define LBMACROFUNCTIONAL

#include "petsc.h"

class LBMacroFunctional {

public:
  virtual void evaluate(PetscInt,const PetscScalar,void*,void*,PetscScalar*) const = 0;
  virtual void adjointCollideSource(PetscInt,const PetscScalar,
                                    void* adjoint, void* fdist, void* obstacle) const = 0;
  virtual void adjointSensitivitySource(PetscInt,const PetscScalar,
                                        void* sens, void* fdist, void* obstacle) const {}
};

#endif
