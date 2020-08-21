#ifndef ADJOINTSYMMETRY2D
#define ADJOINTSYMMETRY2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"

template <class Lattice>
class AdjointSymmetrySouth2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointSymmetrySouth2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointSymmetrySouth2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

#endif
