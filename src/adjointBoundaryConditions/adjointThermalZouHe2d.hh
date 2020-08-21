#ifndef ADJOINTTHERMALZOUHE2D
#define ADJOINTTHERMALZOUHE2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"
#include "latticeBoltzmann/lattices2d.hh"

/************ Adjoint Zou He standard boundaries ************/

/******
       Base templates
******/

template <class IsoLattice, class ThermalLattice>
class AdjointThermalZouHeWest2d : public AdjointLBBoundaryLoop {};

/******
       D2Q4 implementation
******/

template <class IsoLattice>
class AdjointThermalZouHeWest2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalZouHeWest2d(Box2d _box)
    : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tadj;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[0] -= tadj[3];
      tadj[1] -= tadj[3];
      tadj[2] -= tadj[3];
      tadj[3] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

/******
       D2Q5 implementation
******/

template <class IsoLattice>
class AdjointThermalZouHeWest2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalZouHeWest2d(Box2d _box)
    : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tadj;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[0] -= tadj[4];
      tadj[1] -= tadj[4];
      tadj[2] -= tadj[4];
      tadj[3] -= tadj[4];
      tadj[4] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

/************ Adjoint Zou He Neumann boundaries ************/

/******
       Base templates
******/

template <class IsoLattice, class ThermalLattice>
class AdjointThermalZouHeNeumannEast2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalZouHeNeumannWest2d : public AdjointLBBoundaryLoop {};

/******
       D2Q4 implementation
******/

template <class IsoLattice>
class AdjointThermalZouHeNeumannEast2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalZouHeNeumannEast2d(Box2d _box)
    : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt dd;
    PetscScalar* tadj;
    PetscScalar* tadj1;
    PetscScalar* tadj2;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj1 = adj[jj][ii-1] + IsoLattice::numDOF;
      tadj2 = adj[jj][ii-2] + IsoLattice::numDOF;
      for (dd = 0; dd < D2Q4::numDOF; ++dd){
        tadj1[dd] += (4./3.)*tadj[1];
        tadj2[dd] -= (1./3.)*tadj[1];
      }
      tadj[0] -= tadj[1];
      tadj[2] -= tadj[1];
      tadj[3] -= tadj[1];
      tadj[1] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalZouHeNeumannWest2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalZouHeNeumannWest2d(Box2d _box)
    : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt dd;
    PetscScalar* tadj;
    PetscScalar* tadj1;
    PetscScalar* tadj2;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj1 = adj[jj][ii+1] + IsoLattice::numDOF;
      tadj2 = adj[jj][ii+2] + IsoLattice::numDOF;
      for (dd = 0; dd < D2Q4::numDOF; ++dd){
        tadj1[dd] += (4./3.)*tadj[3];
        tadj2[dd] -= (1./3.)*tadj[3];
      }
      tadj[0] -= tadj[3];
      tadj[1] -= tadj[3];
      tadj[2] -= tadj[3];
      tadj[3] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

/******
       D2Q5 implementation
******/

template <class IsoLattice>
class AdjointThermalZouHeNeumannEast2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalZouHeNeumannEast2d(Box2d _box)
    : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt dd;
    PetscScalar* tadj;
    PetscScalar* tadj1;
    PetscScalar* tadj2;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj1 = adj[jj][ii-1] + IsoLattice::numDOF;
      tadj2 = adj[jj][ii-2] + IsoLattice::numDOF;
      for (dd = 0; dd < D2Q5::numDOF; ++dd){
        tadj1[dd] += (4./3.)*tadj[2];
        tadj2[dd] -= (1./3.)*tadj[2];
      }
      tadj[0] -= tadj[2];
      tadj[1] -= tadj[2];
      tadj[3] -= tadj[2];
      tadj[4] -= tadj[2];
      tadj[2] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalZouHeNeumannWest2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalZouHeNeumannWest2d(Box2d _box)
    : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt dd;
    PetscScalar* tadj;
    PetscScalar* tadj1;
    PetscScalar* tadj2;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj1 = adj[jj][ii+1] + IsoLattice::numDOF;
      tadj2 = adj[jj][ii+2] + IsoLattice::numDOF;
      for (dd = 0; dd < D2Q5::numDOF; ++dd){
        tadj1[dd] += (4./3.)*tadj[4];
        tadj2[dd] -= (1./3.)*tadj[4];
      }
      tadj[0] -= tadj[4];
      tadj[1] -= tadj[4];
      tadj[2] -= tadj[4];
      tadj[3] -= tadj[4];
      tadj[4] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

#endif
