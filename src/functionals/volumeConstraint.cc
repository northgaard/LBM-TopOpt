#include "volumeConstraint.hh"

void VolumeConstraint::computeFunctional(Vec domainVec, PetscScalar* value) const
{
  PetscScalar domainSum;
  VecSum(domainVec,&domainSum);

  switch (type){
  case ConstrainFluid:
    *value = -volume + (1./vecSize)*domainSum;
    break;
  case ConstrainMaterial:
    *value = (1. - volume) - (1./vecSize)*domainSum;
    break;
  }
}

void VolumeConstraint::computeSensitivities(Vec domainSensitivities) const
{
  PetscScalar diffValue = 0.;

  switch (type){
  case ConstrainFluid:
    diffValue = 1./vecSize;
    break;
  case ConstrainMaterial:
    diffValue = -1./vecSize;
    break;
  }
  VecSet(domainSensitivities,diffValue);
}

PetscErrorCode VolumeConstraint::computeFluidFraction(Vec domain, PetscScalar* value)
{
  PetscScalar domainSum;
  PetscInt domainSize;
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecSum(domain,&domainSum); CHKERRQ(ierr);
  ierr = VecGetSize(domain,&domainSize); CHKERRQ(ierr);
  *value = domainSum / domainSize; CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VolumeConstraint::computeMaterialFraction(Vec domain, PetscScalar* value)
{
  Vec dummy;
  PetscErrorCode ierr;
  PetscScalar sum;
  PetscInt domainSize;
  PetscFunctionBeginUser;
  ierr = VecGetSize(domain,&domainSize);
  ierr = VecDuplicate(domain,&dummy); CHKERRQ(ierr);
  ierr = VecSet(dummy,1.); CHKERRQ(ierr);
  ierr = VecAXPY(dummy,-1.,domain); CHKERRQ(ierr);
  ierr = VecSum(dummy,&sum); CHKERRQ(ierr);
  *value = sum / domainSize;
  PetscFunctionReturn(0);
}
