#ifndef ADJOINTBOUNCEBACK2D
#define ADJOINTBOUNCEBACK2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"

template <class Lattice>
class AdjointOnGridBouncebackNorth2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackNorth2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackNorth2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackSouth2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackSouth2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackSouth2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackWest2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackWest2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackEast2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackEast2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackNorthWest2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackNorthWest2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackNorthWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackNorthEast2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackNorthEast2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackNorthEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackSouthWest2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackSouthWest2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackSouthWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointOnGridBouncebackSouthEast2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointOnGridBouncebackSouthEast2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointOnGridBouncebackSouthEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

#endif
