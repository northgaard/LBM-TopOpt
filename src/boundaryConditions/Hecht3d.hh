#ifndef HECHT3D
#define HECHT3D

#include "petsc.h"
#include "boundaryConditions/LBBoundaryLoop.hh"
#include "core/geometry3d.hh"
#include "core/macros.hh"
#include "latticeBoltzmann/lattices3d.hh"
#include "boundaryConditions/boundaryDescriptor3d.hh"
#include <type_traits>

/******
       Velocity boundaries -- Base templates
******/

template <class Lattice, class Input>
class StandardHechtVelocityWest3d : public LBBoundaryLoop {};

template <class Lattice, class Input>
class IncompressibleHechtVelocityWest3d : public LBBoundaryLoop {};

/******
       Velocity boundaries -- D3Q19 implementation
******/

template <class Input>
class StandardHechtVelocityWest3d<D3Q19,Input> : public LBBoundaryLoop {

public:
  StandardHechtVelocityWest3d(Box3d _box, Input _vel)
    : LBBoundaryLoop("StandardHechtVelocity3d","WestFace"), boundingBox(_box),
      velFunc(_vel){}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,uy,uz,rho;
    PetscScalar Nyx,Nzx;
    PetscScalar* fd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        fd = fdist[kk][jj][ii];
        velFunc(timestep,ii,jj,kk,ux,uy,uz);
        rho = (1./(1. - ux))*(fd[0]+fd[1]+fd[2]+fd[4]+fd[5]+fd[10]+fd[11]+fd[13]+fd[14]
                              + 2.*(fd[3]+fd[6]+fd[8]+fd[16]+fd[18]));
        Nyx = 0.5*(fd[5]+fd[11]+fd[13]-(fd[2]+fd[4]+fd[14])) - (1./3.)*rho*uy;
        Nzx = 0.5*(fd[10]+fd[13]+fd[14]-(fd[1]+fd[4]+fd[5])) - (1./3.)*rho*uz;
        fd[7] = fd[16] + (rho/6.)*(ux - uz) + Nzx;
        fd[9] = fd[18] + (rho/6.)*(ux - uy) + Nyx;
        fd[12] = fd[3] + (rho/3.)*ux;
        fd[15] = fd[6] + (rho/6.)*(ux + uz) - Nzx;
        fd[17] = fd[8] + (rho/6.)*(ux + uy) - Nyx;
      }
    }
  }
private:
  Box3d boundingBox;
  Input velFunc;
};

template <class Input>
class IncompressibleHechtVelocityWest3d<D3Q19,Input> : public LBBoundaryLoop {

public:
  IncompressibleHechtVelocityWest3d(Box3d _box, Input _vel)
    : LBBoundaryLoop("IncompressibleHechtVelocity3d","WestFace"), boundingBox(_box),
      velFunc(_vel){}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,uy,uz;
    PetscScalar Nyx,Nzx;
    PetscScalar* fd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        fd = fdist[kk][jj][ii];
        velFunc(timestep,ii,jj,kk,ux,uy,uz);
        Nyx = 0.5*(fd[5]+fd[11]+fd[13]-(fd[2]+fd[4]+fd[14])) - (1./3.)*uy;
        Nzx = 0.5*(fd[10]+fd[13]+fd[14]-(fd[1]+fd[4]+fd[5])) - (1./3.)*uz;
        fd[7] = fd[16] + (1./6.)*(ux - uz) + Nzx;
        fd[9] = fd[18] + (1./6.)*(ux - uy) + Nyx;
        fd[12] = fd[3] + (1./3.)*ux;
        fd[15] = fd[6] + (1./6.)*(ux + uz) - Nzx;
        fd[17] = fd[8] + (1./6.)*(ux + uy) - Nyx;
      }
    }
  }
private:
  Box3d boundingBox;
  Input velFunc;
};

/******
       Velocity boundaries -- Boundary descriptor
******/

template <class CollisionOperator, class Lattice, class Input>
class HechtVelocityBoundary3d : public BoundaryDescriptor3d {

public:
  HechtVelocityBoundary3d(Input _vel)
    : BoundaryDescriptor3d("HechtVelocity3d"), velFunc(_vel){}
  PetscErrorCode getWestBoundary(const Box3d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<typename CollisionOperator::Equilibrium,
        StandardFeq3d<Lattice>>::value)
      {
        *loop = new (std::nothrow) StandardHechtVelocityWest3d<Lattice,Input>
          (boundaryLocation,velFunc); CHKNEWPTR(*loop);
      } else {
      *loop = new (std::nothrow) IncompressibleHechtVelocityWest3d<Lattice,Input>
        (boundaryLocation,velFunc); CHKNEWPTR(*loop);
    }
    PetscFunctionReturn(0);
  }
private:
  Input velFunc;
};

/******
       Pressure boundaries -- Base templates
******/

template <class Lattice, class Input>
class StandardHechtPressureEast3d : public LBBoundaryLoop {};

template <class Lattice, class Input>
class IncompressibleHechtPressureEast3d : public LBBoundaryLoop {};

/******
       Pressure boundaries -- D3Q19 implementation
******/

template <class Input>
class StandardHechtPressureEast3d<D3Q19,Input> : public LBBoundaryLoop {

public:
  StandardHechtPressureEast3d(Box3d _box, Input _pres)
    : LBBoundaryLoop("StandardHechtPressure3d","EastFace"), boundingBox(_box),
      presFunc(_pres){}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,rho0;
    PetscScalar Nyx,Nzx;
    PetscScalar* fd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        fd = fdist[kk][jj][ii];
        presFunc(timestep,ii,jj,kk,rho0);
        ux = -1. + (1./rho0)*(fd[0]+fd[1]+fd[2]+fd[4]+fd[5]+fd[10]+fd[11]+fd[13]+fd[14]
                              + 2.*(fd[7]+fd[9]+fd[12]+fd[15]+fd[17]));
        Nyx = 0.5*(fd[5]+fd[11]+fd[13]-(fd[2]+fd[4]+fd[14]));
        Nzx = 0.5*(fd[10]+fd[13]+fd[14]-(fd[1]+fd[4]+fd[5]));
        fd[3] = fd[12] - (1./3.)*rho0*ux;
        fd[6] = fd[15] - (rho0/6.)*ux + Nzx;
        fd[8] = fd[17] - (rho0/6.)*ux + Nyx;
        fd[16] = fd[7] - (rho0/6.)*ux - Nzx;
        fd[18] = fd[9] - (rho0/6.)*ux - Nyx;
      }
    }
  }
private:
  Box3d boundingBox;
  Input presFunc;
};

template <class Input>
class IncompressibleHechtPressureEast3d<D3Q19,Input> : public LBBoundaryLoop {

public:
  IncompressibleHechtPressureEast3d(Box3d _box, Input _pres)
    : LBBoundaryLoop("IncompressibleHechtPressure3d","EastFace"), boundingBox(_box),
      presFunc(_pres){}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,rho0;
    PetscScalar Nyx,Nzx;
    PetscScalar* fd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        fd = fdist[kk][jj][ii];
        presFunc(timestep,ii,jj,kk,rho0);
        ux = fd[0]+fd[1]+fd[2]+fd[4]+fd[5]+fd[10]+fd[11]+fd[13]+fd[14]
          + 2.*(fd[7]+fd[9]+fd[12]+fd[15]+fd[17]) - rho0;
        Nyx = 0.5*(fd[5]+fd[11]+fd[13]-(fd[2]+fd[4]+fd[14]));
        Nzx = 0.5*(fd[10]+fd[13]+fd[14]-(fd[1]+fd[4]+fd[5]));
        fd[3] = fd[12] - (1./3.)*ux;
        fd[6] = fd[15] - (1./6.)*ux + Nzx;
        fd[8] = fd[17] - (1./6.)*ux + Nyx;
        fd[16] = fd[7] - (1./6.)*ux - Nzx;
        fd[18] = fd[9] - (1./6.)*ux - Nyx;
      }
    }
  }
private:
  Box3d boundingBox;
  Input presFunc;
};

/******
       Pressure boundaries -- Boundary descriptor
******/

template <class CollisionOperator, class Lattice, class Input>
class HechtPressureBoundary3d : public BoundaryDescriptor3d {

public:
  HechtPressureBoundary3d(Input _pr)
    : BoundaryDescriptor3d("HechtPressure3d"), presFunc(_pr){}
  PetscErrorCode getEastBoundary(const Box3d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<typename CollisionOperator::Equilibrium,
        StandardFeq3d<Lattice>>::value)
      {
        *loop = new (std::nothrow) StandardHechtPressureEast3d<Lattice,Input>
          (boundaryLocation,presFunc); CHKNEWPTR(*loop);
      } else {
      *loop = new (std::nothrow) IncompressibleHechtPressureEast3d<Lattice,Input>
        (boundaryLocation,presFunc); CHKNEWPTR(*loop);
    }
    PetscFunctionReturn(0);
  }
private:
  Input presFunc;
};

#endif
