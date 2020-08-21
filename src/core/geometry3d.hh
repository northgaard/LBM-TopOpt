#ifndef GEOMETRY3D
#define GEOMETRY3D

#include "petsc.h"
#include "core/globalDefinitions.hh"

struct Box3d {

  BoxRange xRange,yRange,zRange;
  Box3d() : xRange(), yRange(), zRange() {}
  Box3d(PetscInt _x0, PetscInt _x1, PetscInt _y0, PetscInt _y1,
	PetscInt _z0, PetscInt _z1) :
    xRange(_x0,_x1), yRange(_y0,_y1), zRange(_z0,_z1) {}

  PetscInt getNx() const {
    return xRange.getEndId() - xRange.getBeginId() + 1;
  }

  PetscInt getNy() const {
    return yRange.getEndId() - yRange.getBeginId() + 1;
  }

  PetscInt getNz() const {
    return zRange.getEndId() - zRange.getBeginId() + 1;
  }

  void print() const {
    printf("x(%d,%d), y(%d,%d), z(%d,%d)\n",
	   xRange.getBeginId(),xRange.getEndId(),
	   yRange.getBeginId(),yRange.getEndId(),
	   zRange.getBeginId(),zRange.getEndId());
  }

};

struct Box3dBoundaries {
  // Faces
  Box3d north;
  Box3d south;
  Box3d west;
  Box3d east;
  Box3d front;
  Box3d back;
  // Edges
  Box3d northWestEdge;
  Box3d northEastEdge;
  Box3d northFrontEdge;
  Box3d northBackEdge;
  Box3d southWestEdge;
  Box3d southEastEdge;
  Box3d southFrontEdge;
  Box3d southBackEdge;
  Box3d westFrontEdge;
  Box3d westBackEdge;
  Box3d eastFrontEdge;
  Box3d eastBackEdge;
  // Corners
  Box3d northFaceNorthWestCorner;
  Box3d northFaceNorthEastCorner;
  Box3d northFaceSouthWestCorner;
  Box3d northFaceSouthEastCorner;
  Box3d southFaceNorthWestCorner;
  Box3d southFaceNorthEastCorner;
  Box3d southFaceSouthWestCorner;
  Box3d southFaceSouthEastCorner;
};

PetscBool boxIntersection(const Box3d&, const Box3d&, Box3d&);
PetscBool doBoxesIntersect(const Box3d&, const Box3d&);
Box3dBoundaries getBoxBoundaries(const Box3d&);

#endif
