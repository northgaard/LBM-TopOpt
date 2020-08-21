#include "boundaryDescriptor2d.hh"

#undef __FUNCT__
#define __FUNCT__ "getNorthBoundary"
PetscErrorCode BoundaryDescriptor2d::getNorthBoundary(const Box2d,
						      BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "North",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointNorthBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointNorthBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	  "North",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getSouthBoundary"
PetscErrorCode BoundaryDescriptor2d::getSouthBoundary(const Box2d,
						      BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	   "South",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointSouthBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointSouthBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	   "South",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getWestBoundary"
PetscErrorCode BoundaryDescriptor2d::getWestBoundary(const Box2d,
						     BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "West",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointWestBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointWestBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	  "West",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getEastBoundary"
PetscErrorCode BoundaryDescriptor2d::getEastBoundary(const Box2d,
						     BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "East",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointEastBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointEastBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	  "East",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getNorthWestBoundary"
PetscErrorCode BoundaryDescriptor2d::getNorthWestBoundary(const Box2d,
							  BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "NorthWest",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointNorthWestBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointNorthWestBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	  "NorthWest",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getNorthEastBoundary"
PetscErrorCode BoundaryDescriptor2d::getNorthEastBoundary(const Box2d,
							  BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "NorthEast",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointNorthEastBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointNorthEastBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "NorthEast",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getSouthWestBoundary"
PetscErrorCode BoundaryDescriptor2d::getSouthWestBoundary(const Box2d,
							  BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "SouthWest",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointSouthWestBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointSouthWestBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	  "SouthWest",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getSouthEastBoundary"
PetscErrorCode BoundaryDescriptor2d::getSouthEastBoundary(const Box2d,
							  BoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the boundary type %s is not implemented.\n",
	  "SouthEast",boundaryName.c_str());
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getAdjointSouthEastBoundary"
PetscErrorCode BoundaryDescriptor2d::
getAdjointSouthEastBoundary(const Box2d, AdjointBoundaryFunction&) const
{

  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
	  "%s orientation of the adjoint boundary type %s is not implemented.\n",
	  "SouthEast",boundaryName.c_str());
  PetscFunctionReturn(0);

}

/* Refactored implementation starts here */

PetscErrorCode NewBoundaryDescriptor2d::getNorthBoundary(const Box2d,
                                                         LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "North",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getSouthBoundary(const Box2d,
                                                         LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "South",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getWestBoundary(const Box2d,
                                                        LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "West",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getEastBoundary(const Box2d,
                                                        LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "East",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getNorthWestBoundary(const Box2d,
                                                             LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "NorthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getNorthEastBoundary(const Box2d,
                                                             LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "NorthEast",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getSouthWestBoundary(const Box2d,
                                                             LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "SouthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode NewBoundaryDescriptor2d::getSouthEastBoundary(const Box2d,
                                                             LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s orientation of the boundary type %s is not implemented.\n",
           "SouthEast",boundaryName);
  PetscFunctionReturn(0);
}
