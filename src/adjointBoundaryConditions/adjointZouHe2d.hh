#ifndef ADJOINTZOUHE2D
#define ADJOINTZOUHE2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"

/* Velocity boundaries */

template <class Lattice>
class AdjointIncompressibleZouHeVelocityWest2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHeVelocityWest2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHeVelocityWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointIncompressibleZouHeVelocityNorth2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHeVelocityNorth2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHeVelocityNorth2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointIncompressibleZouHeVelocityEast2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHeVelocityEast2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHeVelocityEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice, class Input>
class AdjointStandardZouHeVelocityWest2d : public AdjointLBBoundaryLoop {};

template <class Input>
class AdjointStandardZouHeVelocityWest2d<D2Q9,Input> : public AdjointLBBoundaryLoop {

public:
  AdjointStandardZouHeVelocityWest2d(Box2d _box, Input _vf)
    : velFunc(_vf), boundingBox(_box){}
  void execute(PetscInt timestep, void* adj_v, void*) const override
  {
    PetscScalar*** adj = (PetscScalar***) adj_v;
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,uy;
    PetscScalar adjRho;

    for (auto jj : boundingBox.yRange){
      velFunc(timestep,ii,jj,ux,uy);
      adjRho = ((1./6)*(ux / (1. - ux)) - 0.5*(uy / (1. - ux)))*adj[jj][ii][1]
        + ((1./6)*(ux / (1. - ux)) + 0.5*(uy / (1. - ux)))*adj[jj][ii][7]
        + (2./3.)*(ux / (1. - ux))*adj[jj][ii][8];

      adj[jj][ii][0] += adjRho;
      adj[jj][ii][2] += 0.5*(adj[jj][ii][7] - adj[jj][ii][1]) + adjRho;
      adj[jj][ii][3] += adj[jj][ii][7] + 2.*adjRho;
      adj[jj][ii][4] += adj[jj][ii][8] + 2.*adjRho;
      adj[jj][ii][5] += adj[jj][ii][1] + 2.*adjRho;
      adj[jj][ii][6] += 0.5*(adj[jj][ii][1] - adj[jj][ii][7]) + adjRho;

      adj[jj][ii][1] = 0.;
      adj[jj][ii][7] = 0.;
      adj[jj][ii][8] = 0.;
    }
  }
private:
  Input velFunc;
  Box2d boundingBox;
};

/* Pressure boundaries */

template <class Lattice>
class AdjointIncompressibleZouHePressureEast2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHePressureEast2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHePressureEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointIncompressibleZouHePressureWest2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHePressureWest2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHePressureWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

/* Neumann boundaries */

template <class Lattice>
class AdjointIncompressibleZouHeNeumannEast2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHeNeumannEast2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHeNeumannEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

template <class Lattice>
class AdjointIncompressibleZouHeNeumannWest2d : public AdjointLBBoundaryLoop {};

template <>
class AdjointIncompressibleZouHeNeumannWest2d<D2Q9> : public AdjointLBBoundaryLoop {

public:
  AdjointIncompressibleZouHeNeumannWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt,void*,void*) const override;
private:
  Box2d boundingBox;
};

#endif
