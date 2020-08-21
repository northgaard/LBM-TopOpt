#ifndef VOLUMECONSTRAINT
#define VOLUMECONSTRAINT

#include "petsc.h"
#include "functionals/materialFunctional.hh"

enum VolumeConstraintType { ConstrainFluid, ConstrainMaterial };

class VolumeConstraint : public MaterialFunctional {

public:
  VolumeConstraint(PetscScalar _vol, VolumeConstraintType _t, Vec designDomain)
    : volume(_vol), type(_t)
  {
    VecGetSize(designDomain,&vecSize);
  }
  void computeFunctional(Vec,PetscScalar*) const override;
  void computeSensitivities(Vec) const override;
  void setNewConstraintValue(PetscScalar _v){ volume = _v; }
  static PetscErrorCode computeFluidFraction(Vec,PetscScalar*);
  static PetscErrorCode computeMaterialFraction(Vec,PetscScalar*);
private:
  PetscScalar volume;
  PetscInt vecSize;
  VolumeConstraintType type;
};

#endif
