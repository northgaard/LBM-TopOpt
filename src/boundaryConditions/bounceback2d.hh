#ifndef BOUNCEBACK2D
#define BOUNCEBACK2D

#include "petsc.h"
#include "LBBoundaryLoop.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include <utility>
#include <memory>

template <class Lattice>
class OnGridBouncebackNorth2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackNorth2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorth2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","North"),
                                        boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackSouth2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackSouth2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouth2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","South"),
                                        boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackWest2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackWest2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackWest2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","West"),
                                       boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackEast2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackEast2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackEast2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","East"),
                                       boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackNorthWest2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackNorthWest2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthWest2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","NorthWest"),
                                            boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackNorthEast2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackNorthEast2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthEast2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","NorthEast"),
                                            boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackSouthWest2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackSouthWest2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthWest2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","SouthWest"),
                                            boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class OnGridBouncebackSouthEast2d : public LBBoundaryLoop {};

template <>
class OnGridBouncebackSouthEast2d<D2Q9> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthEast2d(Box2d _box) : LBBoundaryLoop("OnGridBounceback2d","SouthEast"),
                                            boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

/*
  Descriptor
*/

template <class Lattice>
class OnGridBounceback2d : public NewBoundaryDescriptor2d {};

template <>
class OnGridBounceback2d<D2Q9> : public NewBoundaryDescriptor2d {

public:
  OnGridBounceback2d() : NewBoundaryDescriptor2d("OnGridBounceback2d"){}
  PetscErrorCode getNorthBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getSouthBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getWestBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getEastBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getNorthWestBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getNorthEastBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getSouthWestBoundary(const Box2d,LBBoundaryLoop**) const override;
  PetscErrorCode getSouthEastBoundary(const Box2d,LBBoundaryLoop**) const override;
};

#endif
