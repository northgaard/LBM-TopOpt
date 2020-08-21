#ifndef BOUNDARYDESCRIPTOR2D
#define BOUNDARYDESCRIPTOR2D

#include "petsc.h"
#include "LBBoundaryLoop.hh"
#include "core/geometry2d.hh"
#include <functional>
#include <string>

class BoundaryDescriptor2d {

protected:

  using BoundaryFunction = std::function<void(PetscInt,void*)>;
  using AdjointBoundaryFunction = std::function<void(PetscInt,void*,void*)>;

public:

  virtual PetscErrorCode getNorthBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointNorthBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getSouthBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointSouthBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getWestBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointWestBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getEastBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointEastBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getNorthWestBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointNorthWestBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getNorthEastBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointNorthEastBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getSouthWestBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointSouthWestBoundary(const Box2d, AdjointBoundaryFunction&) const;
  virtual PetscErrorCode getSouthEastBoundary(const Box2d,BoundaryFunction&) const;
  virtual PetscErrorCode
  getAdjointSouthEastBoundary(const Box2d, AdjointBoundaryFunction&) const;

  virtual ~BoundaryDescriptor2d(){}

protected:

  BoundaryDescriptor2d(const std::string _bn) : boundaryName(_bn) {}
  std::string boundaryName;

};

class NewBoundaryDescriptor2d {

public:

  virtual PetscErrorCode getNorthBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getWestBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getEastBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthWestBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getNorthEastBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthWestBoundary(const Box2d,LBBoundaryLoop**) const;
  virtual PetscErrorCode getSouthEastBoundary(const Box2d,LBBoundaryLoop**) const;

  virtual ~NewBoundaryDescriptor2d(){}
  char* boundaryName;
protected:
  NewBoundaryDescriptor2d(const char _bn[])
  {
    boundaryName = (char*) _bn;
  }
};

#endif
