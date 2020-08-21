#ifndef INCOMPRESSIBLEMRT2D
#define INCOMPRESSIBLEMRT2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/codiheader.hh"
#include "dynamics/isothermalMacros2d.hh"
#include "dynamics/equilibriumDistributions2d.hh"
#include "adjointDynamics/reverseIncompressibleMRT2d.hh"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
class IncompressibleMRT2d {};

template <class Lattice, class Real = PetscScalar>
class IncompressibleMRTForcing2d {};

template <class Real>
class IncompressibleMRT2d<D2Q9,Real> {

  friend class IncompressibleMRTForcing2d<D2Q9,Real>;
public:
  using LatticeType = D2Q9;
  using RealType = Real;
  using Equilibrium = IncompressibleFeq2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  IncompressibleMRT2d(IncompressibleFlowParameters par){
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  IncompressibleMRT2d(Real _om) : omega(_om){}
  IncompressibleMRT2d() : omega(1.){}
  IncompressibleMRT2d<D2Q9,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return IncompressibleMRT2d<D2Q9,CodiReverseType<Real>>(this->omega);
  }
  ReverseIncompressibleMRT2d<D2Q9,Real> getSourceAdjoint() const
  {
    return ReverseIncompressibleMRT2d<D2Q9,Real>(this->omega);
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         PETSC_RESTRICT Real* mac) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqr = mac[1]*mac[1] + mac[2]*mac[2];
    Real moments[D2Q9::numDOF - 3]; // First three moments are the standard macros
    // Moment transform
    moments[0] = -4.*fdistIn[0] + 2.*fdistIn[1] - fdistIn[2] + 2.*fdistIn[3] -
      fdistIn[4] + 2.*fdistIn[5] - fdistIn[6] + 2.*fdistIn[7] - fdistIn[8];
    moments[1] = 4.*fdistIn[0] + fdistIn[1] - 2.*fdistIn[2] + fdistIn[3] -
      2.*fdistIn[4] + fdistIn[5] - 2.*fdistIn[6] + fdistIn[7] - 2.*fdistIn[8];
    moments[2] = fdistIn[1] - fdistIn[3] + 2.*fdistIn[4] - fdistIn[5] +
      fdistIn[7] - 2.*fdistIn[8];
    moments[3] = -fdistIn[1] + 2.*fdistIn[2] - fdistIn[3] + fdistIn[5] -
      2.*fdistIn[6] + fdistIn[7];
    moments[4] = -fdistIn[2] + fdistIn[4] - fdistIn[6] + fdistIn[8];
    moments[5] = -fdistIn[1] + fdistIn[3] - fdistIn[5] + fdistIn[7];
    // Collision in moment space
    moments[0] = s1*(moments[0] + 2.*mac[0] - 3.*uSqr);
    moments[1] = s2*(moments[1] - mac[0] + 3.*uSqr);
    moments[2] = s46*(moments[2] + mac[1]);
    moments[3] = s46*(moments[3] + mac[2]);
    moments[4] = omega*(moments[4] - mac[1]*mac[1] + mac[2]*mac[2]);
    moments[5] = omega*(moments[5] - mac[1]*mac[2]);
    // Map back to distributions space
    fdistOut[0] = fdistIn[0] + oo9*moments[0] - oo9*moments[1];
    fdistOut[1] = fdistIn[1] - oo18*moments[0] - oo36*moments[1] - oo12*moments[2]
      + oo12*moments[3] + 0.25*moments[5];
    fdistOut[2] = fdistIn[2] + oo36*moments[0] + oo18*moments[1] -
      oo6*moments[3] + 0.25*moments[4];
    fdistOut[3] = fdistIn[3] - oo18*moments[0] - oo36*moments[1] + oo12*moments[2]
      + oo12*moments[3] - 0.25*moments[5];
    fdistOut[4] = fdistIn[4] + oo36*moments[0] + oo18*moments[1] - oo6*moments[2]
      - 0.25*moments[4];
    fdistOut[5] = fdistIn[5] - oo18*moments[0] - oo36*moments[1] + oo12*moments[2]
      - oo12*moments[3] + 0.25*moments[5];
    fdistOut[6] = fdistIn[6] + oo36*moments[0] + oo18*moments[1] + oo6*moments[3]
      + 0.25*moments[4];
    fdistOut[7] = fdistIn[7] - oo18*moments[0] - oo36*moments[1] - oo12*moments[2]
      - oo12*moments[3] - 0.25*moments[5];
    fdistOut[8] = fdistIn[8] + oo36*moments[0] + oo18*moments[1] + oo6*moments[2]
      - 0.25*moments[4];
  }
private:
  static constexpr PetscScalar s1 = 1.4;
  static constexpr PetscScalar s2 = 1.4;
  static constexpr PetscScalar s46 = 1.2;
  Real omega;
  // Fractions
  static constexpr PetscScalar oo6 = 1./6.;
  static constexpr PetscScalar oo9 = 1./9.;
  static constexpr PetscScalar oo12 = 1./12.;
  static constexpr PetscScalar oo18 = 1./18.;
  static constexpr PetscScalar oo36 = 1./36.;
};

template <class Real>
class IncompressibleMRTForcing2d<D2Q9,Real> {

  friend class ReverseIncompressibleMRTForcing2d<D2Q9,Real>;
public:
  using LatticeType = D2Q9;
  using RealType = Real;
  using Equilibrium = IncompressibleFeq2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  IncompressibleMRTForcing2d(IncompressibleFlowParameters par, Real _fx, Real _fy)
    : baseOp(par), omega(baseOp.omega), Fx(_fx), Fy(_fy)
  {}
  IncompressibleMRTForcing2d(Real _om, Real _fx, Real _fy) : baseOp(_om), omega(baseOp.omega),
                                                             Fx(_fx), Fy(_fy){}
  IncompressibleMRTForcing2d() : baseOp(1.), omega(1.), Fx(0.), Fy(0.){}
  IncompressibleMRTForcing2d<D2Q9,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return IncompressibleMRTForcing2d<D2Q9,CodiReverseType<Real>>(this->omega,this->Fx,
                                                           this->Fy);
  }
  ReverseIncompressibleMRTForcing2d<D2Q9,Real> getSourceAdjoint() const
  {
    return ReverseIncompressibleMRTForcing2d<D2Q9,Real>(omega,Fx,Fy);
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         PETSC_RESTRICT Real* mac) const
  {
    baseOp(fdistIn,fdistOut,mac);
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
  static constexpr PetscScalar s1 = IncompressibleMRT2d<D2Q9,Real>::s1;
  static constexpr PetscScalar s2 = IncompressibleMRT2d<D2Q9,Real>::s2;
  static constexpr PetscScalar s46 = IncompressibleMRT2d<D2Q9,Real>::s46;
  IncompressibleMRT2d<D2Q9,Real> baseOp;
  Real omega;
  Real Fx, Fy;
  // Fractions
  static constexpr PetscScalar oo3 = 1./3.;
  static constexpr PetscScalar oo6 = 1./6.;
  static constexpr PetscScalar oo9 = 1./9.;
  static constexpr PetscScalar oo12 = 1./12.;
  static constexpr PetscScalar oo18 = 1./18.;
  static constexpr PetscScalar oo24 = 1./24.;
  static constexpr PetscScalar oo36 = 1./36.;
};

#endif
