#ifndef REVERSEINCOMPRESSIBLEMRT2D
#define REVERSEINCOMPRESSIBLEMRT2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "adjointDynamics/reverseIsothermalMacros2d.hh"

template <class Lattice, class Real = PetscScalar>
class ReverseIncompressibleMRT2d {};

template <class Lattice, class Real = PetscScalar>
class ReverseIncompressibleMRTForcing2d {};

template <class Real>
class ReverseIncompressibleMRT2d<D2Q9,Real> {

public:
  using RevMacros = ReverseIncompressibleMacros2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  using RealType = Real;

  ReverseIncompressibleMRT2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  ReverseIncompressibleMRT2d(Real _om) : omega(_om){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqrRev;
    Real momentsRev[D2Q9::numDOF - 3];
    // Reverse map back to distributions space
    momentsRev[0] = oo36*fdistOutRev[8] - oo18*fdistOutRev[7] + oo36*fdistOutRev[6]
      - oo18*fdistOutRev[5] + oo36*fdistOutRev[4] - oo18*fdistOutRev[3]
      + oo36*fdistOutRev[2] - oo18*fdistOutRev[1] + oo9*fdistOutRev[0];
    momentsRev[1] = oo18*fdistOutRev[8] - oo36*fdistOutRev[7] + oo18*fdistOutRev[6]
      - oo36*fdistOutRev[5] + oo18*fdistOutRev[4] - oo36*fdistOutRev[3]
      + oo18*fdistOutRev[2] - oo36*fdistOutRev[1] - oo9*fdistOutRev[0];
    momentsRev[2] = oo6*fdistOutRev[8] - oo12*fdistOutRev[7] + oo12*fdistOutRev[5]
      - oo6*fdistOutRev[4] + oo12*fdistOutRev[3] - oo12*fdistOutRev[1];
    momentsRev[3] = -oo12*fdistOutRev[7] + oo6*fdistOutRev[6] - oo12*fdistOutRev[5]
      + oo12*fdistOutRev[3] - oo6*fdistOutRev[2] + oo12*fdistOutRev[1];
    momentsRev[4] = -0.25*fdistOutRev[8] + 0.25*fdistOutRev[6] - 0.25*fdistOutRev[4]
      + 0.25*fdistOutRev[2];
    momentsRev[5] = -0.25*fdistOutRev[7] + 0.25*fdistOutRev[5] - 0.25*fdistOutRev[3]
      + 0.25*fdistOutRev[1];
    // Reverse collision in moment space
    momentsRev[5] *= omega;
    momentsRev[4] *= omega;
    momentsRev[3] *= s46;
    momentsRev[2] *= s46;
    momentsRev[1] *= s2;
    momentsRev[0] *= s1;
    macRev[0] += 2.*momentsRev[0] - momentsRev[1];
    macRev[1] += momentsRev[2] - 2.*mac[1]*momentsRev[4] - mac[2]*momentsRev[5];
    macRev[2] += momentsRev[3] + 2.*mac[2]*momentsRev[4] - mac[1]*momentsRev[5];
    uSqrRev = 3.*momentsRev[1] - 3.*momentsRev[0];
    macRev[1] += 2*mac[1]*uSqrRev;
    macRev[2] += 2*mac[2]*uSqrRev;
    // Reverse moment transform
    fdistInRev[0] += fdistOutRev[0] - 4.*momentsRev[0] + 4.*momentsRev[1];
    fdistInRev[1] += fdistOutRev[1] + 2.*momentsRev[0] + momentsRev[1] + momentsRev[2]
      - momentsRev[3] - momentsRev[5];
    fdistInRev[2] += fdistOutRev[2] - momentsRev[0] - 2.*momentsRev[1] + 2.*momentsRev[3]
      - momentsRev[4];
    fdistInRev[3] += fdistOutRev[3] + 2.*momentsRev[0] + momentsRev[1] - momentsRev[2]
      - momentsRev[3] + momentsRev[5];
    fdistInRev[4] += fdistOutRev[4] - momentsRev[0] - 2.*momentsRev[1] + 2.*momentsRev[2]
      + momentsRev[4];
    fdistInRev[5] += fdistOutRev[5] + 2.*momentsRev[0] + momentsRev[1] - momentsRev[2]
      + momentsRev[3] - momentsRev[5];
    fdistInRev[6] += fdistOutRev[6] - momentsRev[0] - 2.*momentsRev[1] - 2.*momentsRev[3]
      - momentsRev[4];
    fdistInRev[7] += fdistOutRev[7] + 2.*momentsRev[0] + momentsRev[1] + momentsRev[2]
      + momentsRev[3] + momentsRev[5];
    fdistInRev[8] += fdistOutRev[8] - momentsRev[0] - 2.*momentsRev[1] - 2.*momentsRev[2]
      + momentsRev[4];
    RevMacros::compute(fdistInRev,macRev[0],macRev[1],macRev[2]);
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
class ReverseIncompressibleMRTForcing2d<D2Q9,Real> {

public:
  using RevMacros = ReverseIncompressibleMacros2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  using RealType = Real;

  ReverseIncompressibleMRTForcing2d(Real _om, Real _fx, Real _fy)
    : reverseBaseOp(_om), omega(_om), Fx(_fx), Fy(_fy){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Real tempb, tempb0, tempb1, tempb2, tempb3, tempb4, tempb5, tempb6, tempb7;
    macRev[1] = macRev[1] + oo12*Fx*(s1-s2*2.-omega*3.+8.)*fdistOutRev[8];
    macRev[2] = macRev[2] + (s1-s2*2.+omega*3.-4.)*oo12*Fy*fdistOutRev[8];
    tempb = oo24*Fx*fdistOutRev[7];
    tempb0 = -(oo6*Fy*fdistOutRev[7]);
    macRev[1] = macRev[1] + (omega*0.75-1.5)*tempb0 + (s2*2.-s1*4.+4.)*tempb;
    macRev[2] = macRev[2] + (s1-s2*0.5-1.)*tempb0 + (6.-omega*3.)*tempb;
    macRev[2] = macRev[2] + oo12*Fy*(s1-s2*2.-omega*3.+8.)*fdistOutRev[6];
    macRev[1] = macRev[1] + (s1-s2*2.+omega*3.-4.)*oo12*Fx*fdistOutRev[6];
    tempb1 = oo24*Fx*fdistOutRev[5];
    tempb2 = -(oo6*Fy*fdistOutRev[5]);
    macRev[1] = macRev[1] + (1.5-omega*0.75)*tempb2 + (s2*2.-s1*4.+4.)*tempb1;
    macRev[2] = macRev[2] + (s1-s2*0.5-1.)*tempb2 + (omega*3.-6.)*tempb1;
    macRev[1] = macRev[1] + oo12*Fx*(s1-s2*2.-omega*3.+8.)*fdistOutRev[4];
    macRev[2] = macRev[2] + (s1-s2*2.+omega*3.-4.)*oo12*Fy*fdistOutRev[4];
    tempb3 = oo24*Fx*fdistOutRev[3];
    tempb4 = -(oo6*Fy*fdistOutRev[3]);
    macRev[1] = macRev[1] + (omega*0.75-1.5)*tempb4 + (s2*2.-s1*4.+4.)*tempb3;
    macRev[2] = macRev[2] + (s1-s2*0.5-1.)*tempb4 + (6.-omega*3.)*tempb3;
    macRev[2] = macRev[2] + oo12*Fy*(s1-s2*2.-omega*3.+8.)*fdistOutRev[2];
    macRev[1] = macRev[1] + (s1-s2*2.+omega*3.-4.)*oo12*Fx*fdistOutRev[2];
    tempb5 = oo24*Fx*fdistOutRev[1];
    tempb6 = -(oo6*Fy*fdistOutRev[1]);
    macRev[1] = macRev[1] + (1.5-omega*0.75)*tempb6 + (s2*2.-s1*4.+4.)*tempb5;
    macRev[2] = macRev[2] + (s1-s2*0.5-1.)*tempb6 + (omega*3.-6.)*tempb5;
    tempb7 = oo3*(s1+s2-4.)*fdistOutRev[0];
    macRev[1] = macRev[1] + Fx*tempb7;
    macRev[2] = macRev[2] + Fy*tempb7;
    reverseBaseOp(fdistIn,fdistOut,mac,fdistInRev,fdistOutRev,macRev);
  }
private:
  static constexpr PetscScalar s1 = 1.4;
  static constexpr PetscScalar s2 = 1.4;
  static constexpr PetscScalar s46 = 1.2;
  ReverseIncompressibleMRT2d<D2Q9,Real> reverseBaseOp;
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
