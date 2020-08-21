#include "collisionDescriptorBases2d.hh"

#undef __FUNCT__
#define __FUNCT__ "getAdjointStreamingLoop"
PetscErrorCode ObstacleCollisionDescriptorBase2d::
getAdjointStreamingLoop(Box2d,Box2d,BaseDynamicsFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "Adjoint dynamics are not implemented for the specified collision descriptor.\n"
	  );
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointCollisionLoop"
PetscErrorCode ObstacleCollisionDescriptorBase2d::
getAdjointCollisionLoop(Box2d,AdjointObstacleDynamicsFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "Adjoint dynamics are not implemented for the specified collision descriptor.\n"
	  );
  PetscFunctionReturn(0);

}
