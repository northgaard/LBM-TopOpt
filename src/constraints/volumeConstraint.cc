#include "volumeConstraint.hh"

PetscErrorCode VolumeConstraint::computeConstraint(Vec domainVec, PetscScalar& value)
{

  PetscErrorCode ierr;
  PetscScalar domainSum;
  PetscInt vecSize;

  ierr = VecSum(domainVec,&domainSum); CHKERRQ(ierr);
  ierr = VecGetSize(domainVec,&vecSize); CHKERRQ(ierr);

  switch (type) {

  case ConstrainFluid:
    value = -volume + (1./vecSize)*domainSum;
    break;
  case ConstrainMaterial:
    value = (1. - volume) - (1./vecSize)*domainSum;
    break;

  }

  return 0;

}

PetscErrorCode VolumeConstraint::computeSensitivities(Vec domainSensitivities)
{

  PetscErrorCode ierr;
  PetscScalar diffValue = 0.;
  PetscInt vecSize;

  ierr = VecGetSize(domainSensitivities,&vecSize); CHKERRQ(ierr);

  switch (type) {

  case ConstrainFluid:
    diffValue = 1./vecSize;
    break;
  case ConstrainMaterial:
    diffValue = -1./vecSize;
    break;

  }

  ierr = VecSet(domainSensitivities,diffValue); CHKERRQ(ierr);

  return 0;

}
