#ifndef BOUNDARYDESCRIPTOR3D
#define BOUNDARYDESCRIPTOR3D

#include "petsc.h"
#include "LBBoundaryLoop.hh"
#include "core/geometry3d.hh"

class BoundaryDescriptor3d {

public:
  /* Faces */
  virtual PetscErrorCode getNorthBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getWestBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getEastBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getFrontBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getBackBoundary(const Box3d,LBBoundaryLoop**) const;
  /* Edges */
  virtual PetscErrorCode getNorthWestEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthEastEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthFrontEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthBackEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthWestEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthEastEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthFrontEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthBackEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getWestFrontEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getWestBackEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getEastFrontEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getEastBackEdgeBoundary(const Box3d,LBBoundaryLoop**) const;
  /* Corners */
  virtual PetscErrorCode getNorthFaceNorthWestCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthFaceNorthEastCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthFaceSouthWestCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthFaceSouthEastCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthFaceNorthWestCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthFaceNorthEastCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthFaceSouthWestCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthFaceSouthEastCornerBoundary(const Box3d,
                                                             LBBoundaryLoop**) const;

  virtual ~BoundaryDescriptor3d(){}
  char* boundaryName;
protected:
  BoundaryDescriptor3d(const char _bn[])
  {
    boundaryName = (char*) _bn;
  }
};

#endif
