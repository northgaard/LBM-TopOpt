#include "geometry3d.hh"

PetscBool boxIntersection(const Box3d& b1, const Box3d& b2, Box3d& theBox)
{

  PetscInt x0,x1,y0,y1,z0,z1;
  x0 = PetscMax(b1.xRange.getBeginId(),b2.xRange.getBeginId());
  x1 = PetscMin(b1.xRange.getEndId(),b2.xRange.getEndId());
  y0 = PetscMax(b1.yRange.getBeginId(),b2.yRange.getBeginId());
  y1 = PetscMin(b1.yRange.getEndId(),b2.yRange.getEndId());
  z0 = PetscMax(b1.zRange.getBeginId(),b2.zRange.getBeginId());
  z1 = PetscMin(b1.zRange.getEndId(),b2.zRange.getEndId());

  theBox = Box3d(x0,x1,y0,y1,z0,z1);
  PetscBool doesIntersect = (PetscBool) (x1 >= x0 && y1 >= y0 && z1 >= z0);

  return doesIntersect;

}

PetscBool doBoxesIntersect(const Box3d& b1, const Box3d& b2)
{

  PetscInt x0,x1,y0,y1,z0,z1;
  x0 = PetscMax(b1.xRange.getBeginId(),b2.xRange.getBeginId());
  x1 = PetscMin(b1.xRange.getEndId(),b2.xRange.getEndId());
  y0 = PetscMax(b1.yRange.getBeginId(),b2.yRange.getBeginId());
  y1 = PetscMin(b1.yRange.getEndId(),b2.yRange.getEndId());
  z0 = PetscMax(b1.zRange.getBeginId(),b2.xRange.getBeginId());
  z1 = PetscMin(b1.zRange.getEndId(),b2.zRange.getEndId());

  return (PetscBool) (x1 >= x0 && y1 >= y0 && z1 >= z0);

}

Box3dBoundaries getBoxBoundaries(const Box3d& theBox)
{

  PetscInt x0,x1,y0,y1,z0,z1;
  Box3dBoundaries boxB;
  x0 = theBox.xRange.getBeginId();
  y0 = theBox.yRange.getBeginId();
  z0 = theBox.zRange.getBeginId();
  x1 = theBox.xRange.getEndId();
  y1 = theBox.yRange.getEndId();
  z1 = theBox.zRange.getEndId();
  // Faces
  boxB.north = Box3d(x0+1,x1-1,y1,y1,z0+1,z1-1);
  boxB.south = Box3d(x0+1,x1-1,y0,y0,z0+1,z1-1);
  boxB.west = Box3d(x0,x0,y0+1,y1-1,z0+1,z1-1);
  boxB.east = Box3d(x1,x1,y0+1,y1-1,z0+1,z1-1);
  boxB.front = Box3d(x0+1,x1-1,y0+1,y1-1,z1,z1);
  boxB.back = Box3d(x0+1,x1-1,y0+1,y1-1,z0,z0);
  // Edges
  boxB.northWestEdge = Box3d(x0,x0,y1,y1,z0+1,z1-1);
  boxB.northEastEdge = Box3d(x1,x1,y1,y1,z0+1,z1-1);
  boxB.northFrontEdge = Box3d(x0+1,x1-1,y1,y1,z1,z1);
  boxB.northBackEdge = Box3d(x0+1,x1-1,y1,y1,z0,z0);
  boxB.southWestEdge = Box3d(x0,x0,y0,y0,z0+1,z1-1);
  boxB.southEastEdge = Box3d(x1,x1,y0,y0,z0+1,z1-1);
  boxB.southFrontEdge = Box3d(x0+1,x1-1,y0,y0,z1,z1);
  boxB.southBackEdge = Box3d(x0+1,x1-1,y0,y0,z0,z0);
  boxB.westFrontEdge = Box3d(x0,x0,y0+1,y1-1,z1,z1);
  boxB.westBackEdge = Box3d(x0,x0,y0+1,y1-1,z0,z0);
  boxB.eastFrontEdge = Box3d(x1,x1,y0+1,y1-1,z1,z1);
  boxB.eastBackEdge = Box3d(x1,x1,y0+1,y1-1,z0,z0);
  // Corners
  boxB.northFaceNorthWestCorner = Box3d(x0,x0,y1,y1,z0,z0);
  boxB.northFaceNorthEastCorner = Box3d(x1,x1,y1,y1,z0,z0);
  boxB.northFaceSouthWestCorner = Box3d(x0,x0,y1,y1,z1,z1);
  boxB.northFaceSouthEastCorner = Box3d(x1,x1,y1,y1,z1,z1);
  boxB.southFaceNorthWestCorner = Box3d(x0,x0,y0,y0,z0,z0);
  boxB.southFaceNorthEastCorner = Box3d(x1,x1,y0,y0,z0,z0);
  boxB.southFaceSouthWestCorner = Box3d(x0,x0,y0,y0,z1,z1);
  boxB.southFaceSouthEastCorner = Box3d(x1,x1,y0,y0,z1,z1);

  return boxB;
}
