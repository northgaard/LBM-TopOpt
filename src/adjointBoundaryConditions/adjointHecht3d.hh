#ifndef ADJOINTHECHT3D
#define ADJOINTHECHT3D

#include "petsc.h"
#include "core/geometry3d.hh"
#include "latticeBoltzmann/lattices3d.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"

/******
       Adjoint velocity boundaries -- Base templates
******/

template <class Lattice>
class AdjointIncompressibleHechtVelocityWest3d : public AdjointLBBoundaryLoop {};

/******
       Adjoint velocity boundaries -- D3Q19 implementation
******/

template <>
class AdjointIncompressibleHechtVelocityWest3d<D3Q19> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleHechtVelocityWest3d(Box3d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box3d boundingBox;
};

/******
       Adjoint pressure boundaries -- Base templates
******/

template <class Lattice>
class AdjointIncompressibleHechtPressureEast3d<Lattice> : public AdjointLBBoundaryLoop {};

/******
       Adjoint pressure boundaries -- Base templates
******/

template <>
class AdjointIncompressibleHechtPressureEast3d<D3Q19> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleHechtPressureEast3d(Box3d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box3d boundingBox;
};

#endif
