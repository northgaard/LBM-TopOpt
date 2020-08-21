#include "geometry2d.hh"

PetscBool boxIntersection(const Box2d& b1, const Box2d& b2, Box2d& theBox)
{

  PetscInt x0,x1,y0,y1;
  x0 = PetscMax(b1.xRange.getBeginId(),b2.xRange.getBeginId());
  x1 = PetscMin(b1.xRange.getEndId(),b2.xRange.getEndId());
  y0 = PetscMax(b1.yRange.getBeginId(),b2.yRange.getBeginId());
  y1 = PetscMin(b1.yRange.getEndId(),b2.yRange.getEndId());

  theBox = Box2d(x0,x1,y0,y1);
  PetscBool doesIntersect = (PetscBool) (x1 >= x0 && y1 >= y0);

  return doesIntersect;

}

PetscBool doBoxesIntersect(const Box2d& b1, const Box2d& b2)
{

  PetscInt x0,x1,y0,y1;
  x0 = PetscMax(b1.xRange.getBeginId(),b2.xRange.getBeginId());
  x1 = PetscMin(b1.xRange.getEndId(),b2.xRange.getEndId());
  y0 = PetscMax(b1.yRange.getBeginId(),b2.yRange.getBeginId());
  y1 = PetscMin(b1.yRange.getEndId(),b2.yRange.getEndId());

  return (PetscBool) (x1 >= x0 && y1 >= y0);
}

Box2dBoundaries getBoxBoundaries(const Box2d& theBox)
{

  PetscInt x0,x1,y0,y1;
  Box2dBoundaries boxB;
  x0 = theBox.xRange.getBeginId();
  y0 = theBox.yRange.getBeginId();
  x1 = theBox.xRange.getEndId();
  y1 = theBox.yRange.getEndId();
  
  boxB.north = Box2d(x0+1,x1-1,y1,y1);
  boxB.south = Box2d(x0+1,x1-1,y0,y0);
  boxB.west = Box2d(x0,x0,y0+1,y1-1);
  boxB.east = Box2d(x1,x1,y0+1,y1-1);

  boxB.northWest = Box2d(x0,x0,y1,y1);
  boxB.northEast = Box2d(x1,x1,y1,y1);
  boxB.southWest = Box2d(x0,x0,y0,y0);
  boxB.southEast = Box2d(x1,x1,y0,y0);

  return boxB;

}
