#ifndef GEOMETRY2D
#define GEOMETRY2D

#include "petsc.h"
#include "core/globalDefinitions.hh"

struct Box2d {

  BoxRange xRange,yRange;
  Box2d() : xRange(), yRange() {}
  Box2d(PetscInt _x0, PetscInt _x1, PetscInt _y0, PetscInt _y1) :
    xRange(_x0,_x1), yRange(_y0,_y1) {}

  PetscInt getNx() const {
    return xRange.getEndId() - xRange.getBeginId() + 1;
  }

  PetscInt getNy() const {
    return yRange.getEndId() - yRange.getBeginId() + 1;
  }

  void print() const {
    printf("x(%d,%d), y(%d,%d)\n",xRange.getBeginId(),xRange.getEndId(),
	   yRange.getBeginId(),yRange.getEndId());
  }

};

struct Box2dBoundaries {

  // Edges
  Box2d north;
  Box2d south;
  Box2d west;
  Box2d east;
  // Corners
  Box2d northWest;
  Box2d northEast;
  Box2d southWest;
  Box2d southEast;

};

PetscBool boxIntersection(const Box2d&, const Box2d&, Box2d&);
PetscBool doBoxesIntersect(const Box2d&, const Box2d&);
Box2dBoundaries getBoxBoundaries(const Box2d&);

#endif
