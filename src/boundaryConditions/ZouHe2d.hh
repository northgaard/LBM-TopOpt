#ifndef ZOUHE2D
#define ZOUHE2D

#include "petsc.h"
#include "boundaryConditions/LBBoundaryLoop.hh"
#include "core/geometry2d.hh"
#include "core/macros.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "adjointBoundaryConditions/adjointZouHe2d.hh"
#include "dynamics/equilibriumDistributions2d.hh"
#include <type_traits>

/*
  Velocity boundaries
*/

/* Incompressible */

template <class Lattice, class Input>
class IncompressibleZouHeVelocityWest2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleZouHeVelocityWest2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleZouHeVelocityWest2d(Box2d _box, Input _vel)
    : LBBoundaryLoop("IncompressibleZouHeVelocity2d","West"), boundingBox(_box), velFunc(_vel){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,uy;

    for (auto jj : boundingBox.yRange){
      velFunc(timestep,ii,jj,ux,uy);
      fdistOut[jj][ii][1] = 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]) + fdistIn[jj][ii][5]
        + (1./6.)*ux - 0.5*uy;
      fdistOut[jj][ii][7] = 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]) + fdistIn[jj][ii][3]
        + (1./6.)*ux + 0.5*uy;
      fdistOut[jj][ii][8] = fdistIn[jj][ii][4] + (2./3.)*ux;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHeVelocityWest2d<D2Q9>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input velFunc;
};

template <class Lattice, class Input>
class IncompressibleZouHeVelocityNorth2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleZouHeVelocityNorth2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleZouHeVelocityNorth2d(Box2d _box, Input _vel)
    : LBBoundaryLoop("IncompressibleZouHeVelocity2d","North"), boundingBox(_box), velFunc(_vel){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar ux,uy;

    for (auto ii : boundingBox.xRange){
      velFunc(timestep,ii,jj,ux,uy);
      fdistOut[jj][ii][1] = 0.5*(fdistIn[jj][ii][4] - fdistIn[jj][ii][8]) + fdistIn[jj][ii][5]
        - (1./6.)*uy + 0.5*ux;
      fdistOut[jj][ii][2] = fdistIn[jj][ii][6] - (2./3.)*uy;
      fdistOut[jj][ii][3] = 0.5*(fdistIn[jj][ii][8] - fdistIn[jj][ii][4]) + fdistIn[jj][ii][7]
        - (1./6.)*uy - 0.5*ux;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHeVelocityNorth2d<D2Q9>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input velFunc;
};

template <class Lattice, class Input>
class IncompressibleZouHeVelocityEast2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleZouHeVelocityEast2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleZouHeVelocityEast2d(Box2d _box, Input _vel)
    : LBBoundaryLoop("IncompressibleZouHeVelocity2d","East"), boundingBox(_box), velFunc(_vel){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux,uy;

    for (auto jj : boundingBox.yRange){
      velFunc(timestep,ii,jj,ux,uy);
      fdistOut[jj][ii][3] = fdistIn[jj][ii][7]
        + 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]) - (1./6.)*ux - 0.5*uy;
      fdistOut[jj][ii][4] = fdistIn[jj][ii][8] - (2./3.)*ux;
      fdistOut[jj][ii][5] = fdistIn[jj][ii][1]
        + 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]) - (1./6.)*ux + 0.5*uy;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHeVelocityEast2d<D2Q9>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input velFunc;
};

/* Standard */

template <class Lattice, class Input>
class StandardZouHeVelocityWest2d : public LBBoundaryLoop {};

template <class Input>
class StandardZouHeVelocityWest2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  StandardZouHeVelocityWest2d(Box2d _box, Input _vel)
    : LBBoundaryLoop("StandardZouHeVelocity2d","West"), boundingBox(_box), velFunc(_vel){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar rho,ux,uy;

    for (auto jj : boundingBox.yRange){
      velFunc(timestep,ii,jj,ux,uy);
      rho = (fdistIn[jj][ii][0] + fdistIn[jj][ii][2] + fdistIn[jj][ii][6] +
             2.*(fdistIn[jj][ii][3] + fdistIn[jj][ii][4] + fdistIn[jj][ii][5])) /
        (1. - ux);

      fdistOut[jj][ii][1] = fdistIn[jj][ii][5] +
        0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]) +
        (1./6.)*rho*ux - 0.5*rho*uy;
      fdistOut[jj][ii][7] = fdistIn[jj][ii][3] +
        0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]) +
        (1./6.)*rho*ux + 0.5*rho*uy;
      fdistOut[jj][ii][8] = fdistIn[jj][ii][4] + (2./3.)*rho*ux;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointStandardZouHeVelocityWest2d<D2Q9,Input>
      (boundingBox,velFunc);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input velFunc;
};

template <class CollisionOperator, class Lattice, class Input>
class ZouHeVelocity2d : public NewBoundaryDescriptor2d {};

template <class CollisionOperator, class Input>
class ZouHeVelocity2d<CollisionOperator,D2Q9,Input> : public NewBoundaryDescriptor2d {

public:
  ZouHeVelocity2d(Input _vel)
    : NewBoundaryDescriptor2d("ZouHeVelocity2d"), velFunc(_vel){}
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,StandardFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) StandardZouHeVelocityWest2d<D2Q9,Input>
        (boundaryLocation,velFunc);
      CHKNEWPTR(*loop);
    } else if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHeVelocityWest2d<D2Q9,Input>
        (boundaryLocation,velFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He velocity boundary.\n");
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHeVelocityNorth2d<D2Q9,Input>
        (boundaryLocation,velFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He velocity boundary.\n");
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHeVelocityEast2d<D2Q9,Input>
        (boundaryLocation,velFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He velocity boundary.\n");
    }
    PetscFunctionReturn(0);
  }
private:
  using Equilibrium = typename CollisionOperator::Equilibrium;
  Input velFunc;
};

/*
  Pressure boundaries
*/

// Incompressible
template <class Lattice, class Input>
class IncompressibleZouHePressureEast2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleZouHePressureEast2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleZouHePressureEast2d(Box2d _box, Input _inp)
    : LBBoundaryLoop("IncompressibleZouHePressure2d","East"), boundingBox(_box),
      presFunc(_inp){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux0,rho0;

    for (auto jj : boundingBox.yRange){
      presFunc(timestep,ii,jj,rho0);
      ux0 = -rho0 + (fdistIn[jj][ii][0] + fdistIn[jj][ii][2] + fdistIn[jj][ii][6]
                     + 2.*(fdistIn[jj][ii][1] + fdistIn[jj][ii][7] + fdistIn[jj][ii][8]));

      fdistOut[jj][ii][4] = fdistIn[jj][ii][8] - (2./3.)*ux0;
      fdistOut[jj][ii][3] = fdistIn[jj][ii][7] - (1./6.)*ux0
        + 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]);
      fdistOut[jj][ii][5] = fdistIn[jj][ii][1] - (1./6.)*ux0
        + 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]);
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHePressureEast2d<D2Q9>(boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input presFunc;
};

template <class Lattice, class Input>
class IncompressibleZouHePressureWest2d : public LBBoundaryLoop {};

template <class Input>
class IncompressibleZouHePressureWest2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  IncompressibleZouHePressureWest2d(Box2d _box, Input _inp)
    : LBBoundaryLoop("IncompressibleZouHePressure2d","West"), boundingBox(_box),
      presFunc(_inp){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux0,rho0;

    for (auto jj : boundingBox.yRange){
      presFunc(timestep,ii,jj,rho0);
      ux0 = rho0 - (fdistIn[jj][ii][0] + fdistIn[jj][ii][2] + fdistIn[jj][ii][6]
                    + 2.*(fdistIn[jj][ii][3] + fdistIn[jj][ii][4] + fdistIn[jj][ii][5]));

      fdistOut[jj][ii][1] = fdistIn[jj][ii][5] + (1./6.)*ux0
        + 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]);
      fdistOut[jj][ii][7] = fdistIn[jj][ii][3] + (1./6.)*ux0
        + 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]);
      fdistOut[jj][ii][8] = fdistIn[jj][ii][4] + (2./3.)*ux0;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHePressureWest2d<D2Q9>(boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input presFunc;
};

// Standard
template <class Lattice, class Input>
class StandardZouHePressureEast2d : public LBBoundaryLoop {};

template <class Input>
class StandardZouHePressureEast2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  StandardZouHePressureEast2d(Box2d _box, Input _inp)
    : LBBoundaryLoop("StandardZouHePressure2d","East"), boundingBox(_box),
      presFunc(_inp){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux0,rho0;

    for (auto jj : boundingBox.yRange){
      presFunc(timestep,ii,jj,rho0);
      ux0 = -1. + (fdistIn[jj][ii][0] + fdistIn[jj][ii][2] + fdistIn[jj][ii][6]
                   + 2.*(fdistIn[jj][ii][1] + fdistIn[jj][ii][7] + fdistIn[jj][ii][8])) / rho0;

      fdistOut[jj][ii][4] = fdistIn[jj][ii][8] - (2./3.)*rho0*ux0;
      fdistOut[jj][ii][3] = fdistIn[jj][ii][7] - (1./6.)*rho0*ux0
        + 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]);
      fdistOut[jj][ii][5] = fdistIn[jj][ii][1] - (1./6.)*rho0*ux0
        + 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]);
    }
  }
private:
  Box2d boundingBox;
  Input presFunc;
};

template <class Lattice, class Input>
class StandardZouHePressureWest2d : public LBBoundaryLoop {};

template <class Input>
class StandardZouHePressureWest2d<D2Q9,Input> : public LBBoundaryLoop {

public:
  StandardZouHePressureWest2d(Box2d _box, Input _inp)
    : LBBoundaryLoop("StandardZouHePressure2d","West"), boundingBox(_box),
      presFunc(_inp){}
  void execute(PetscInt timestep,void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux0,rho0;

    for (auto jj : boundingBox.yRange){
      presFunc(timestep,ii,jj,rho0);
      ux0 = 1. - (fdistIn[jj][ii][0] + fdistIn[jj][ii][2] + fdistIn[jj][ii][6]
                  + 2.*(fdistIn[jj][ii][3] + fdistIn[jj][ii][4] + fdistIn[jj][ii][5])) / rho0;

      fdistOut[jj][ii][1] = fdistIn[jj][ii][5] + (1./6.)*rho0*ux0
        + 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]);
      fdistOut[jj][ii][7] = fdistIn[jj][ii][3] + (1./6.)*rho0*ux0
        + 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]);
      fdistOut[jj][ii][8] = fdistIn[jj][ii][4] + (2./3.)*rho0*ux0;
    }
  }
private:
  Box2d boundingBox;
  Input presFunc;
};

template <class Lattice, class Equilibrium, class Input>
class ZouHePressure2d : public NewBoundaryDescriptor2d {};

template <class Equilibrium,class Input>
class ZouHePressure2d<D2Q9,Equilibrium,Input> : public NewBoundaryDescriptor2d {

public:
  ZouHePressure2d(Input _inp)
    : NewBoundaryDescriptor2d("ZouHePressure2d"), presFunc(_inp){}
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHePressureEast2d<D2Q9,Input>
        (boundaryLocation,presFunc);
      CHKNEWPTR(*loop);
    } else if (std::is_same<Equilibrium,StandardFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) StandardZouHePressureEast2d<D2Q9,Input>
        (boundaryLocation,presFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He pressure boundary.\n");
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHePressureWest2d<D2Q9,Input>
        (boundaryLocation,presFunc);
      CHKNEWPTR(*loop);
    } else if (std::is_same<Equilibrium,StandardFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) StandardZouHePressureWest2d<D2Q9,Input>
        (boundaryLocation,presFunc);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He pressure boundary.\n");
    }
    PetscFunctionReturn(0);
  }
private:
  Input presFunc;
};

/*
  Neumann boundaries
*/

template <class Lattice>
class IncompressibleZouHeNeumannEast2d : public LBBoundaryLoop {};

template <class Lattice>
class IncompressibleZouHeNeumannWest2d : public LBBoundaryLoop {};

template <>
class IncompressibleZouHeNeumannEast2d<D2Q9> : public LBBoundaryLoop {

public:
  IncompressibleZouHeNeumannEast2d(Box2d _box)
    : LBBoundaryLoop("IncompressibleZouHeNeumann2d","East"), boundingBox(_box){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux, uxN1 = 0., uxN2 = 0.;
    const PetscScalar* fdistN;

    for (auto jj : boundingBox.yRange){
      fdistN = fdistIn[jj][ii-1];
      Macros::computeUx(fdistN,uxN1);
      fdistN = fdistIn[jj][ii-2];
      Macros::computeUx(fdistN,uxN2);
      ux = (4.*uxN1 - uxN2)/3.;
      fdistOut[jj][ii][3] = fdistIn[jj][ii][7] + 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2])
        - (1./6.)*ux;
      fdistOut[jj][ii][4] = fdistIn[jj][ii][8] - (2./3.)*ux;
      fdistOut[jj][ii][5] = fdistIn[jj][ii][1] + 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6])
        - (1./6.)*ux;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHeNeumannEast2d<D2Q9>(boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  using Macros = IncompressibleMacros2d<D2Q9>;
  Box2d boundingBox;
};

template <>
class IncompressibleZouHeNeumannWest2d<D2Q9> : public LBBoundaryLoop {

public:
  IncompressibleZouHeNeumannWest2d(Box2d _box)
    : LBBoundaryLoop("IncompressibleZouHeNeumann2d","West"), boundingBox(_box){}
  void execute(PetscInt timestep, void* fdistIn_v, void* fdistOut_v) const override
  {
    auto fdistIn = static_cast<const PetscScalar***>(fdistIn_v);
    auto fdistOut = static_cast<PetscScalar***>(fdistOut_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar ux, uxN1 = 0., uxN2 = 0.;
    const PetscScalar* fdistN;

    for (auto jj : boundingBox.yRange){
      fdistN = fdistIn[jj][ii+1];
      Macros::computeUx(fdistN,uxN1);
      fdistN = fdistIn[jj][ii+2];
      Macros::computeUx(fdistN,uxN2);
      ux = (4.*uxN1 - uxN2)/3.;
      fdistOut[jj][ii][1] = 0.5*(fdistIn[jj][ii][6] - fdistIn[jj][ii][2]) + fdistIn[jj][ii][5]
        + (1./6.)*ux;
      fdistOut[jj][ii][7] = 0.5*(fdistIn[jj][ii][2] - fdistIn[jj][ii][6]) + fdistIn[jj][ii][3]
        + (1./6.)*ux;
      fdistOut[jj][ii][8] = fdistIn[jj][ii][4] + (2./3.)*ux;
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointIncompressibleZouHeNeumannWest2d<D2Q9>(boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  using Macros = IncompressibleMacros2d<D2Q9>;
  Box2d boundingBox;
};

template <class CollisionOperator, class Lattice>
class ZouHeNeumann2d : public NewBoundaryDescriptor2d {};

template <class CollisionOperator>
class ZouHeNeumann2d<CollisionOperator,D2Q9> : public NewBoundaryDescriptor2d {

public:
  ZouHeNeumann2d() : NewBoundaryDescriptor2d("ZouHeNeumann2d"){}
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHeNeumannEast2d<D2Q9>
        (boundaryLocation);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He Neumann boundary.\n");
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    if (std::is_same<Equilibrium,IncompressibleFeq2d<D2Q9>>::value){
      *loop = new (std::nothrow) IncompressibleZouHeNeumannWest2d<D2Q9>
        (boundaryLocation);
      CHKNEWPTR(*loop);
    } else {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unrecognized or incompatible equilibrium distribution for Zou He Neumann boundary.\n");
    }
    PetscFunctionReturn(0);
  }
private:
  using Equilibrium = typename CollisionOperator::Equilibrium;
};

#endif
