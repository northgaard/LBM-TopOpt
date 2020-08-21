#ifndef DARCYMRT2D
#define DARCYMRT2D

#include "petsc.h"
#include "dynamics/incompressibleMRT2d.hh"

template <class Lattice, class InterpolationFunction, class Real = PetscScalar>
class DarcyMRT2d {};

template <class InterpolationFunction, class Real>
class DarcyMRT2d<D2Q9,InterpolationFunction,Real> {

public:
  using LatticeType = D2Q9;
  using RealType = Real;
  using Equilibrium = IncompressibleFeq2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  static constexpr PetscInt numAdditionalFields = 1;

  DarcyMRT2d(Real _om, Real _fx, Real _fy, Real _am, InterpolationFunction _inter)
    : baseOp(_om), interp(_inter), omega(_om), cfx(_fx), cfy(_fy), amax(_am) {}
  DarcyMRT2d() : baseOp(), interp(), omega(1.), cfx(0.), cfy(0.), amax(0.) {}
  DarcyMRT2d<D2Q9,CodiReverseType<Real>> getCodiAdjoint() const
  {
    auto codiInter = interp.getCodiAdjoint();
    return DarcyMRT2d<D2Q9,decltype(codiInter),CodiReverseType<Real>>(this->omega,
                                                                      this->cfx,
                                                                      this->cfy,
                                                                      codiInter);
  }
  void getSourceAdjoint() const {}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  PETSC_RESTRICT Real* mac, PETSC_RESTRICT Real* obst)
  {
    Real imult, Fx, Fy;
    baseOp(fdistIn,fdistOut,mac);
    imult = -amax*interp(*obst)*mac[0];
    Fx = imult*mac[1] + cfx;
    Fy = imult*mac[2] + cfy;
    // Add forcing
    fdistOut[0] += oo3*(s1 + s2 - 4.)*(Fx*mac[1] + Fy*mac[2]);
    fdistOut[1] += oo24*((-4.*s1 + 2.*s2 + 4.)*mac[1] + 2. + (3.*omega - 6.)*mac[2])*Fx
      - oo6*Fy*((-0.75*omega + 1.5)*mac[1] + 0.5 + (s1 - 0.5*s2 - 1.)*mac[2]);
    fdistOut[2] += oo12*(-4. + (s1 - 2.*s2 - 3.*omega + 8.)*mac[2])*Fy
      + oo12*Fx*mac[1]*(s1 - 2.*s2 + 3.*omega - 4.);
    fdistOut[3] += oo24*((-4.*s1 + 2.*s2 + 4.)*mac[1] - 2. + (-3.*omega + 6.)*mac[2])*Fx
      - oo6*Fy*((0.75*omega - 1.5)*mac[1] + 0.5 + (s1 - 0.5*s2 - 1.)*mac[2]);
    fdistOut[4] += oo12*(-4. + (s1 - 2.*s2 - 3.*omega + 8.)*mac[1])*Fx
      + oo12*Fy*mac[2]*(s1 - 2.*s2 + 3.*omega - 4.);
    fdistOut[5] += oo24*((-4.*s1 + 2.*s2 + 4.)*mac[1] - 2. + (3.*omega - 6.)*mac[2])*Fx
      - oo6*Fy*((-0.75*omega + 1.5)*mac[1] - 0.5 + (s1 - 0.5*s2 - 1.)*mac[2]);
    fdistOut[6] += oo12*(4. + (s1 - 2.*s2 - 3.*omega + 8.)*mac[2])*Fy
      + oo12*Fx*mac[1]*(s1 - 2.*s2 + 3.*omega - 4.);
    fdistOut[7] += oo24*((-4.*s1 + 2.*s2 + 4.)*mac[1] + 2. + (-3.*omega + 6.)*mac[2])*Fx
      - oo6*Fy*((0.75*omega - 1.5)*mac[1] - 0.5 + (s1 - 0.5*s2 - 1.)*mac[2]);
    fdistOut[8] += oo12*(4. + (s1 - 2.*s2 - 3.*omega + 8.)*mac[1])*Fx
      + oo12*Fy*mac[2]*(s1 - 2.*s2 + 3.*omega - 4.);
  }
private:
  static constexpr PetscScalar s1 = 1.4;
  static constexpr PetscScalar s2 = 1.4;
  static constexpr PetscScalar s46 = 1.2;
  IncompressibleMRT2d<D2Q9,Real> baseOp;
  InterpolationFunction interp;
  Real omega;
  const Real cfx, cfy, amax;
  static constexpr PetscScalar oo3 = 1./3.;
  static constexpr PetscScalar oo6 = 1./6.;
  static constexpr PetscScalar oo9 = 1./9.;
  static constexpr PetscScalar oo12 = 1./12.;
  static constexpr PetscScalar oo18 = 1./18.;
  static constexpr PetscScalar oo24 = 1./24.;
  static constexpr PetscScalar oo36 = 1./36.;
};

#endif
