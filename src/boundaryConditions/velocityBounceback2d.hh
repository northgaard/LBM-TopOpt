#ifndef VELOCITYBOUNCEBACK2D
#define VELOCITYBOUNCEBACK2D

#include "petsc.h"
#include "LBBoundaryLoop.hh"
#include "core/geometry2d.hh"
#include "core/macros.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"

template <class Lattice, class Input>
class IncompressibleBouncebackVelocityWest2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleBouncebackVelocityWest2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleBouncebackVelocityWest2d(Box2d _box, Input _vel)
    : LBBoundaryLoop("BouncebackVelocity2d","West"), boundingBox(_box), velFunc(_vel){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,uy;

    for (auto jj : boundingBox.yRange){
      velFunc(timestep,ii,jj,ux,uy);
      fdistOut[jj][ii][1] = fdistIn[jj][ii][1]
        + 2.*D2Q9::csSqInv*D2Q9::weights[1]*(D2Q9::ex[1]*ux + D2Q9::ey[1]*uy);
      fdistOut[jj][ii][7] = fdistIn[jj][ii][7]
        + 2.*D2Q9::csSqInv*D2Q9::weights[7]*(D2Q9::ex[7]*ux + D2Q9::ey[7]*uy);
      fdistOut[jj][ii][8] = fdistIn[jj][ii][8]
        + 2.*D2Q9::csSqInv*D2Q9::weights[8]*(D2Q9::ex[8]*ux + D2Q9::ey[8]*uy);
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    *adjLoop = nullptr;
    return 0;
  }
private:
  Box2d boundingBox;
  Input velFunc;
};

template <class Lattice, class Input>
class IncompressibleBouncebackVelocityNorth2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleBouncebackVelocityNorth2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleBouncebackVelocityNorth2d(Box2d _box, Input _vel)
    : LBBoundaryLoop("BouncebackVelocity2d","North"), boundingBox(_box), velFunc(_vel){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar ux,uy;

    for (auto ii : boundingBox.xRange){
      velFunc(timestep,ii,jj,ux,uy);
      fdistOut[jj][ii][1] = fdistIn[jj][ii][1]
        + 2.*D2Q9::csSqInv*D2Q9::weights[1]*(D2Q9::ex[1]*ux + D2Q9::ey[1]*uy);
      fdistOut[jj][ii][2] = fdistIn[jj][ii][2]
        + 2.*D2Q9::csSqInv*D2Q9::weights[2]*(D2Q9::ex[2]*ux + D2Q9::ey[2]*uy);
      fdistOut[jj][ii][3] = fdistIn[jj][ii][3]
        + 2.*D2Q9::csSqInv*D2Q9::weights[3]*(D2Q9::ex[3]*ux + D2Q9::ey[3]*uy);
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    *adjLoop = nullptr;
    return 0;
  }
private:
  Box2d boundingBox;
  Input velFunc;
};

template <class CollisionOperator, class Lattice, class Input>
class BouncebackVelocity2d : public NewBoundaryDescriptor2d {};

template <class CollisionOperator, class Input>
class BouncebackVelocity2d<CollisionOperator,D2Q9,Input> : public NewBoundaryDescriptor2d {

public:
  BouncebackVelocity2d(Input _vel)
    : NewBoundaryDescriptor2d("BouncebackVelocity2d"), velFunc(_vel){}
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleBouncebackVelocityWest2d<D2Q9,Input>
        (boundaryLocation,velFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for bounceback velocity boundary.\n");
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleBouncebackVelocityNorth2d<D2Q9,Input>
        (boundaryLocation,velFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for bounceback velocity boundary.\n");
    }
    PetscFunctionReturn(0);
  }
private:
  using Equilibrium = typename CollisionOperator::Equilibrium;
  Input velFunc;
};

#endif
