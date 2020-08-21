#ifndef VOLUMECONSTRAINT
#define VOLUMECONSTRAINT

#include "petsc.h"

enum VolumeConstraintType { ConstrainFluid, ConstrainMaterial };

class VolumeConstraint {

public:

  VolumeConstraint(PetscScalar _vol, VolumeConstraintType _t)
    : volume(_vol), type(_t){}
  VolumeConstraint() : volume(1.), type(ConstrainFluid){}

  PetscErrorCode computeConstraint(Vec,PetscScalar&);
  PetscErrorCode computeSensitivities(Vec);
  
private:

  PetscScalar volume;
  VolumeConstraintType type;

};

#endif
