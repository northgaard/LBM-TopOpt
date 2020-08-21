#ifndef SYMMETRY2D
#define SYMMETRY2D

#include "petsc.h"
#include "LBBoundaryLoop.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"

template <class Lattice>
class SymmetrySouth2d : public LBBoundaryLoop {};

template <>
class SymmetrySouth2d<D2Q9> : public LBBoundaryLoop {

public:
  SymmetrySouth2d(Box2d _box) : LBBoundaryLoop("Symmetry2d","South"),
                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class Symmetry2d : public NewBoundaryDescriptor2d {};

template <>
class Symmetry2d<D2Q9> : public NewBoundaryDescriptor2d {

public:
  Symmetry2d() : NewBoundaryDescriptor2d("Symmetry2d"){}
  PetscErrorCode getSouthBoundary(const Box2d,LBBoundaryLoop**) const override;
};

#endif
