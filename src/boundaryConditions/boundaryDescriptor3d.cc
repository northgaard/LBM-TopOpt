#include "boundaryDescriptor3d.hh"

/***
    Faces
***/

PetscErrorCode BoundaryDescriptor3d::getNorthBoundary(const Box3d,
                                                      LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s face orientation of the boundary type %s is not implemented.\n",
           "North",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getSouthBoundary(const Box3d,
                                                      LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s face orientation of the boundary type %s is not implemented.\n",
           "South",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getWestBoundary(const Box3d,
                                                     LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s face orientation of the boundary type %s is not implemented.\n",
           "West",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getEastBoundary(const Box3d,
                                                     LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s face orientation of the boundary type %s is not implemented.\n",
           "East",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getFrontBoundary(const Box3d,
                                                      LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s face orientation of the boundary type %s is not implemented.\n",
           "Front",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getBackBoundary(const Box3d,
                                                     LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s face orientation of the boundary type %s is not implemented.\n",
           "Back",boundaryName);
  PetscFunctionReturn(0);
}

/***
    Edges
***/

PetscErrorCode BoundaryDescriptor3d::getNorthWestEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "NorthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getNorthEastEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "NorthEast",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getNorthFrontEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "NorthFront",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getNorthBackEdgeBoundary(const Box3d,
                                                               LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "NorthBack",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getSouthWestEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "SouthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getSouthEastEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "SouthEast",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getSouthFrontEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "SouthFront",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getSouthBackEdgeBoundary(const Box3d,
                                                               LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "SouthBack",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getWestFrontEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "WestFront",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getWestBackEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "WestBack",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getEastFrontEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "EastFront",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::getEastBackEdgeBoundary(const Box3d,
                                                              LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s edge orientation of the boundary type %s is not implemented.\n",
           "EastBack",boundaryName);
  PetscFunctionReturn(0);
}

/***
    Corners
***/

PetscErrorCode BoundaryDescriptor3d::
getNorthFaceNorthWestCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "NorthFaceNorthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getNorthFaceNorthEastCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "NorthFaceNorthEast",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getNorthFaceSouthWestCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "NorthFaceSouthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getNorthFaceSouthEastCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "NorthFaceSouthEast",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getSouthFaceNorthWestCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "SouthFaceNorthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getSouthFaceNorthEastCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "SouthFaceNorthEast",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getSouthFaceSouthWestCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "SouthFaceSouthWest",boundaryName);
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryDescriptor3d::
getSouthFaceSouthEastCornerBoundary(const Box3d,LBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "%s corner orientation of the boundary type %s is not implemented.\n",
           "SouthFaceSouthEast",boundaryName);
  PetscFunctionReturn(0);
}
