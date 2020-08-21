#ifndef MRTCOLLISION3D
#define MRTCOLLISION3D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "dynamics/isothermalMacros3d.hh"
#include "latticeBoltzmann/lattices3d.hh"

struct MRT3d_Compile_constants {
  static constexpr PetscScalar s1 = 1.19;
  static constexpr PetscScalar s2 = 1.4;
  static constexpr PetscScalar s4 = 1.2;
  static constexpr PetscScalar s10 = 1.4;
  /*
    s16 = 1.98 is the value given in the paper by d'Humières. I have not found this
    to be stable when velocity boundary conditions come into play.
  */
  static constexpr PetscScalar s16 = 1.4;
  static constexpr PetscScalar we = 3.;
  static constexpr PetscScalar wej = -5.5;
  static constexpr PetscScalar wxx = -0.5;
};

template <class Lattice, class Real = PetscScalar,
          class Constants = MRT3d_Compile_constants>
class IncompressibleMRT3d {};

template <class Real, class Constants>
class IncompressibleMRT3d<D3Q19,Real,Constants> {

public:
  using RealType = Real;
  using Equilibrium = IncompressibleFeq3d<D3Q19,Real>;
  using Macros = IncompressibleMacros3d<D3Q19,Real>;

  IncompressibleMRT3d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  IncompressibleMRT3d(Real _om) : omega(_om){}
  IncompressibleMRT3d() : omega(1.){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac)
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2],mac[3]);
    Real moments[D3Q19::numDOF - 4];
    // Moment transform
    moments[0] = -30.*fdistIn[0] - 11.*fdistIn[1] - 11.*fdistIn[2] - 11.*fdistIn[3]
      + 8.*fdistIn[4] + 8.*fdistIn[5] + 8.*fdistIn[6] + 8.*fdistIn[7]
      + 8.*fdistIn[8] + 8.*fdistIn[9] - 11.*fdistIn[10] - 11.*fdistIn[11]
      - 11.*fdistIn[12] + 8.*fdistIn[13] + 8.*fdistIn[14] + 8.*fdistIn[15]
      + 8.*fdistIn[16] + 8.*fdistIn[17] + 8.*fdistIn[18];
    moments[1] = 12.*fdistIn[0] - 4.*fdistIn[1] - 4.*fdistIn[2] - 4.*fdistIn[3]
      + fdistIn[4] + fdistIn[5] + fdistIn[6] + fdistIn[7] + fdistIn[8]
      + fdistIn[9] - 4.*fdistIn[10] - 4.*fdistIn[11] - 4.*fdistIn[12]
      + fdistIn[13] + fdistIn[14] + fdistIn[15] + fdistIn[16] + fdistIn[17]
      + fdistIn[18];
    moments[2] = 4.*fdistIn[3] - fdistIn[6] + fdistIn[7] - fdistIn[8]
      + fdistIn[9] - 4.*fdistIn[12] + fdistIn[15] - fdistIn[16] + fdistIn[17]
      - fdistIn[18];
    moments[3] = 4.*fdistIn[2] - fdistIn[4] + fdistIn[5] - fdistIn[8] - fdistIn[9]
      - 4.*fdistIn[11] + fdistIn[13] - fdistIn[14] + fdistIn[17] + fdistIn[18];
    moments[4] = 4.*fdistIn[1] - fdistIn[4] - fdistIn[5] - fdistIn[6]
      - fdistIn[7] - 4.*fdistIn[10] + fdistIn[13] + fdistIn[14] + fdistIn[15]
      + fdistIn[16];
    moments[5] = -fdistIn[1] - fdistIn[2] + 2.*fdistIn[3] - 2.*fdistIn[4]
      - 2.*fdistIn[5] + fdistIn[6] + fdistIn[7] + fdistIn[8] + fdistIn[9]
      - fdistIn[10] - fdistIn[11] + 2.*fdistIn[12] - 2.*fdistIn[13]
      - 2.*fdistIn[14] + fdistIn[15] + fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[6] = 2.*fdistIn[1] + 2.*fdistIn[2] - 4.*fdistIn[3] - 2.*fdistIn[4]
      - 2.*fdistIn[5] + fdistIn[6] + fdistIn[7] + fdistIn[8] + fdistIn[9]
      + 2.*fdistIn[10] + 2.*fdistIn[11] - 4.*fdistIn[12] - 2.*fdistIn[13]
      - 2.*fdistIn[14] + fdistIn[15] + fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[7] = -fdistIn[1] + fdistIn[2] - fdistIn[6] - fdistIn[7]
      + fdistIn[8] + fdistIn[9] - fdistIn[10] + fdistIn[11] - fdistIn[15]
      - fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[8] = 2.*fdistIn[1] - 2.*fdistIn[2] - fdistIn[6] - fdistIn[7]
      + fdistIn[8] + fdistIn[9] + 2.*fdistIn[10] - 2.*fdistIn[11]
      - fdistIn[15] - fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[9] = fdistIn[8] - fdistIn[9] + fdistIn[17] - fdistIn[18];
    moments[10] = fdistIn[4] - fdistIn[5] + fdistIn[13] - fdistIn[14];
    moments[11] = fdistIn[6] - fdistIn[7] + fdistIn[15] - fdistIn[16];
    moments[12] = fdistIn[6] - fdistIn[7] - fdistIn[8] + fdistIn[9]
      - fdistIn[15] + fdistIn[16] + fdistIn[17] - fdistIn[18];
    moments[13] = -fdistIn[4] + fdistIn[5] + fdistIn[8] + fdistIn[9]
      + fdistIn[13] - fdistIn[14] - fdistIn[17] - fdistIn[18];
    moments[14] = fdistIn[4] + fdistIn[5] - fdistIn[6] - fdistIn[7]
      - fdistIn[13] - fdistIn[14] + fdistIn[15] + fdistIn[16];
    // Collision in moment space
    Real jSqr = mac[1]*mac[1] + mac[2]*mac[2] + mac[3]*mac[3];
    moments[0] = s1*(moments[0] + 11.*mac[0] - 19.*jSqr);
    moments[1] = s2*(moments[1] - we*mac[0] - wej*jSqr);
    moments[2] = s4*(moments[2] + (2./3.)*mac[1]);
    moments[3] = s4*(moments[3] + (2./3.)*mac[2]);
    moments[4] = s4*(moments[4] + (2./3.)*mac[3]);
    Real tempEq = (2.*mac[1]*mac[1] - (mac[2]*mac[2] + mac[3]*mac[3]));
    moments[5] = omega*(moments[5] - tempEq);
    moments[6] = s10*(moments[6] - wxx*tempEq);
    tempEq = (mac[2]*mac[2] - mac[3]*mac[3]);
    moments[7] = omega*(moments[7] - tempEq);
    moments[8] = s10*(moments[8] - wxx*tempEq);
    moments[9] = omega*(moments[9] - mac[1]*mac[2]);
    moments[10] = omega*(moments[10] - mac[2]*mac[3]);
    moments[11] = omega*(moments[11] - mac[1]*mac[3]);
    moments[12] *= s16;
    moments[13] *= s16;
    moments[14] *= s16;
    // Map back to distribution space
    fdistOut[0] = fdistIn[0] + fo399*moments[0] - oo21*moments[1];
    fdistOut[1] = fdistIn[1] + eo2394*moments[0] + oo63*moments[1]
      - 0.1*moments[4] + oo36*moments[5] - oo36*moments[6]
      + oo12*moments[7] - oo12*moments[8];
    fdistOut[2] = fdistIn[2] + eo2394*moments[0] + oo63*moments[1]
      - 0.1*moments[3] + oo36*moments[5] - oo36*moments[6]
      - oo12*moments[7] + oo12*moments[8];
    fdistOut[3] = fdistIn[3] + eo2394*moments[0] + oo63*moments[1]
      - 0.1*moments[2] - oo18*moments[5] + oo18*moments[6];
    fdistOut[4] = fdistIn[4] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[3] + oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] - 0.25*moments[10] + 0.125*moments[13]
      - 0.125*moments[14];
    fdistOut[5] = fdistIn[5] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[3] + oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] + 0.25*moments[10] - 0.125*moments[13]
      - 0.125*moments[14];
    fdistOut[6] = fdistIn[6] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] + oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      - 0.25*moments[11] - 0.125*moments[12] + 0.125*moments[14];
    fdistOut[7] = fdistIn[7] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] + oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      + 0.25*moments[11] + 0.125*moments[12] + 0.125*moments[14];
    fdistOut[8] = fdistIn[8] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] + oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      - 0.25*moments[9] + 0.125*moments[12] - 0.125*moments[13];
    fdistOut[9] = fdistIn[9] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] + oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      + 0.25*moments[9] - 0.125*moments[12] - 0.125*moments[13];
    fdistOut[10] = fdistIn[10] + eo2394*moments[0] + oo63*moments[1]
      + 0.1*moments[4] + oo36*moments[5] - oo36*moments[6]
      + oo12*moments[7] - oo12*moments[8];
    fdistOut[11] = fdistIn[11] + eo2394*moments[0] + oo63*moments[1]
      + 0.1*moments[3] + oo36*moments[5] - oo36*moments[6]
      - oo12*moments[7] + oo12*moments[8];
    fdistOut[12] = fdistIn[12] + eo2394*moments[0] + oo63*moments[1]
      + 0.1*moments[2] - oo18*moments[5] + oo18*moments[6];
    fdistOut[13] = fdistIn[13] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[3] - oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] - 0.25*moments[10] - 0.125*moments[13]
      + 0.125*moments[14];
    fdistOut[14] = fdistIn[14] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[3] - oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] + 0.25*moments[10] + 0.125*moments[13]
      + 0.125*moments[14];
    fdistOut[15] = fdistIn[15] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] - oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      - 0.25*moments[11] + 0.125*moments[12] - 0.125*moments[14];
    fdistOut[16] = fdistIn[16] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] - oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      + 0.25*moments[11] - 0.125*moments[12] - 0.125*moments[14];
    fdistOut[17] = fdistIn[17] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] - oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      - 0.25*moments[9] - 0.125*moments[12] + 0.125*moments[13];
    fdistOut[18] = fdistIn[18] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] - oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      + 0.25*moments[9] + 0.125*moments[12] + 0.125*moments[13];
  }
private:
  Real omega;
  static constexpr PetscScalar s1 = Constants::s1;
  static constexpr PetscScalar s2 = Constants::s2;
  static constexpr PetscScalar s4 = Constants::s4;
  static constexpr PetscScalar s10 = Constants::s10;
  static constexpr PetscScalar s16 = Constants::s16;
  static constexpr PetscScalar we = Constants::we;
  static constexpr PetscScalar wej = Constants::wej;
  static constexpr PetscScalar wxx = Constants::wxx;
  // Fractions
  static constexpr PetscScalar fo399 = 5./399.;
  static constexpr PetscScalar eo2394 = 11./2394.;
  static constexpr PetscScalar fo1197 = 4./1197.;
  static constexpr PetscScalar oo21 = 1./21.;
  static constexpr PetscScalar oo63 = 1./63.;
  static constexpr PetscScalar oo36 = 1./36.;
  static constexpr PetscScalar oo12 = 1./12.;
  static constexpr PetscScalar oo18 = 1./18.;
  static constexpr PetscScalar oo40 = 1./40.;
  static constexpr PetscScalar oo252 = 1./252.;
  static constexpr PetscScalar oo72 = 1./72.;
  static constexpr PetscScalar oo24 = 1./24.;
};

template <class Lattice, class Real = PetscScalar,
          class Constants = MRT3d_Compile_constants>
class StandardMRT3d {};

template <class Real, class Constants>
class StandardMRT3d<D3Q19,Real,Constants> {

public:
  using RealType = Real;
  using Equilibrium = StandardFeq3d<D3Q19,Real>;
  using Macros = StandardMacros3d<D3Q19,Real>;

  StandardMRT3d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  StandardMRT3d(Real _om) : omega(_om){}
  StandardMRT3d() : omega(1.){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac)
  {
    Macros::computeMoments(fdistIn,mac[0],mac[1],mac[2],mac[3]);
    Real moments[D3Q19::numDOF - 4];
    // Moment transform
    moments[0] = -30.*fdistIn[0] - 11.*fdistIn[1] - 11.*fdistIn[2] - 11.*fdistIn[3]
      + 8.*fdistIn[4] + 8.*fdistIn[5] + 8.*fdistIn[6] + 8.*fdistIn[7]
      + 8.*fdistIn[8] + 8.*fdistIn[9] - 11.*fdistIn[10] - 11.*fdistIn[11]
      - 11.*fdistIn[12] + 8.*fdistIn[13] + 8.*fdistIn[14] + 8.*fdistIn[15]
      + 8.*fdistIn[16] + 8.*fdistIn[17] + 8.*fdistIn[18];
    moments[1] = 12.*fdistIn[0] - 4.*fdistIn[1] - 4.*fdistIn[2] - 4.*fdistIn[3]
      + fdistIn[4] + fdistIn[5] + fdistIn[6] + fdistIn[7] + fdistIn[8]
      + fdistIn[9] - 4.*fdistIn[10] - 4.*fdistIn[11] - 4.*fdistIn[12]
      + fdistIn[13] + fdistIn[14] + fdistIn[15] + fdistIn[16] + fdistIn[17]
      + fdistIn[18];
    moments[2] = 4.*fdistIn[3] - fdistIn[6] + fdistIn[7] - fdistIn[8]
      + fdistIn[9] - 4.*fdistIn[12] + fdistIn[15] - fdistIn[16] + fdistIn[17]
      - fdistIn[18];
    moments[3] = 4.*fdistIn[2] - fdistIn[4] + fdistIn[5] - fdistIn[8] - fdistIn[9]
      - 4.*fdistIn[11] + fdistIn[13] - fdistIn[14] + fdistIn[17] + fdistIn[18];
    moments[4] = 4.*fdistIn[1] - fdistIn[4] - fdistIn[5] - fdistIn[6]
      - fdistIn[7] - 4.*fdistIn[10] + fdistIn[13] + fdistIn[14] + fdistIn[15]
      + fdistIn[16];
    moments[5] = -fdistIn[1] - fdistIn[2] + 2.*fdistIn[3] - 2.*fdistIn[4]
      - 2.*fdistIn[5] + fdistIn[6] + fdistIn[7] + fdistIn[8] + fdistIn[9]
      - fdistIn[10] - fdistIn[11] + 2.*fdistIn[12] - 2.*fdistIn[13]
      - 2.*fdistIn[14] + fdistIn[15] + fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[6] = 2.*fdistIn[1] + 2.*fdistIn[2] - 4.*fdistIn[3] - 2.*fdistIn[4]
      - 2.*fdistIn[5] + fdistIn[6] + fdistIn[7] + fdistIn[8] + fdistIn[9]
      + 2.*fdistIn[10] + 2.*fdistIn[11] - 4.*fdistIn[12] - 2.*fdistIn[13]
      - 2.*fdistIn[14] + fdistIn[15] + fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[7] = -fdistIn[1] + fdistIn[2] - fdistIn[6] - fdistIn[7]
      + fdistIn[8] + fdistIn[9] - fdistIn[10] + fdistIn[11] - fdistIn[15]
      - fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[8] = 2.*fdistIn[1] - 2.*fdistIn[2] - fdistIn[6] - fdistIn[7]
      + fdistIn[8] + fdistIn[9] + 2.*fdistIn[10] - 2.*fdistIn[11]
      - fdistIn[15] - fdistIn[16] + fdistIn[17] + fdistIn[18];
    moments[9] = fdistIn[8] - fdistIn[9] + fdistIn[17] - fdistIn[18];
    moments[10] = fdistIn[4] - fdistIn[5] + fdistIn[13] - fdistIn[14];
    moments[11] = fdistIn[6] - fdistIn[7] + fdistIn[15] - fdistIn[16];
    moments[12] = fdistIn[6] - fdistIn[7] - fdistIn[8] + fdistIn[9]
      - fdistIn[15] + fdistIn[16] + fdistIn[17] - fdistIn[18];
    moments[13] = -fdistIn[4] + fdistIn[5] + fdistIn[8] + fdistIn[9]
      + fdistIn[13] - fdistIn[14] - fdistIn[17] - fdistIn[18];
    moments[14] = fdistIn[4] + fdistIn[5] - fdistIn[6] - fdistIn[7]
      - fdistIn[13] - fdistIn[14] + fdistIn[15] + fdistIn[16];
    // Collision in moment space
    Real jSqr = mac[1]*mac[1] + mac[2]*mac[2] + mac[3]*mac[3];
    PetscScalar rhoInv = 1./mac[0];
    moments[0] = s1*(moments[0] + 11.*mac[0] - 19.*rhoInv*jSqr);
    moments[1] = s2*(moments[1] - we*mac[0] - wej*rhoInv*jSqr);
    moments[2] = s4*(moments[2] + (2./3.)*mac[1]);
    moments[3] = s4*(moments[3] + (2./3.)*mac[2]);
    moments[4] = s4*(moments[4] + (2./3.)*mac[3]);
    Real tempEq = rhoInv*(2.*mac[1]*mac[1] - (mac[2]*mac[2] + mac[3]*mac[3]));
    moments[5] = omega*(moments[5] - tempEq);
    moments[6] = s10*(moments[6] - wxx*tempEq);
    tempEq = rhoInv*(mac[2]*mac[2] - mac[3]*mac[3]);
    moments[7] = omega*(moments[7] - tempEq);
    moments[8] = s10*(moments[8] - wxx*tempEq);
    moments[9] = omega*(moments[9] - rhoInv*mac[1]*mac[2]);
    moments[10] = omega*(moments[10] - rhoInv*mac[2]*mac[3]);
    moments[11] = omega*(moments[11] - rhoInv*mac[1]*mac[3]);
    moments[12] *= s16;
    moments[13] *= s16;
    moments[14] *= s16;
    // Map back to distribution space
    fdistOut[0] = fdistIn[0] + fo399*moments[0] - oo21*moments[1];
    fdistOut[1] = fdistIn[1] + eo2394*moments[0] + oo63*moments[1]
      - 0.1*moments[4] + oo36*moments[5] - oo36*moments[6]
      + oo12*moments[7] - oo12*moments[8];
    fdistOut[2] = fdistIn[2] + eo2394*moments[0] + oo63*moments[1]
      - 0.1*moments[3] + oo36*moments[5] - oo36*moments[6]
      - oo12*moments[7] + oo12*moments[8];
    fdistOut[3] = fdistIn[3] + eo2394*moments[0] + oo63*moments[1]
      - 0.1*moments[2] - oo18*moments[5] + oo18*moments[6];
    fdistOut[4] = fdistIn[4] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[3] + oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] - 0.25*moments[10] + 0.125*moments[13]
      - 0.125*moments[14];
    fdistOut[5] = fdistIn[5] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[3] + oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] + 0.25*moments[10] - 0.125*moments[13]
      - 0.125*moments[14];
    fdistOut[6] = fdistIn[6] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] + oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      - 0.25*moments[11] - 0.125*moments[12] + 0.125*moments[14];
    fdistOut[7] = fdistIn[7] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] + oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      + 0.25*moments[11] + 0.125*moments[12] + 0.125*moments[14];
    fdistOut[8] = fdistIn[8] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] + oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      - 0.25*moments[9] + 0.125*moments[12] - 0.125*moments[13];
    fdistOut[9] = fdistIn[9] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] + oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      + 0.25*moments[9] - 0.125*moments[12] - 0.125*moments[13];
    fdistOut[10] = fdistIn[10] + eo2394*moments[0] + oo63*moments[1]
      + 0.1*moments[4] + oo36*moments[5] - oo36*moments[6]
      + oo12*moments[7] - oo12*moments[8];
    fdistOut[11] = fdistIn[11] + eo2394*moments[0] + oo63*moments[1]
      + 0.1*moments[3] + oo36*moments[5] - oo36*moments[6]
      - oo12*moments[7] + oo12*moments[8];
    fdistOut[12] = fdistIn[12] + eo2394*moments[0] + oo63*moments[1]
      + 0.1*moments[2] - oo18*moments[5] + oo18*moments[6];
    fdistOut[13] = fdistIn[13] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[3] - oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] - 0.25*moments[10] - 0.125*moments[13]
      + 0.125*moments[14];
    fdistOut[14] = fdistIn[14] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[3] - oo40*moments[4] + oo18*moments[5]
      + oo36*moments[6] + 0.25*moments[10] + 0.125*moments[13]
      + 0.125*moments[14];
    fdistOut[15] = fdistIn[15] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] - oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      - 0.25*moments[11] + 0.125*moments[12] - 0.125*moments[14];
    fdistOut[16] = fdistIn[16] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] - oo40*moments[4] - oo36*moments[5]
      - oo72*moments[6] + oo12*moments[7] + oo24*moments[8]
      + 0.25*moments[11] - 0.125*moments[12] - 0.125*moments[14];
    fdistOut[17] = fdistIn[17] - fo1197*moments[0] - oo252*moments[1]
      - oo40*moments[2] - oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      - 0.25*moments[9] - 0.125*moments[12] + 0.125*moments[13];
    fdistOut[18] = fdistIn[18] - fo1197*moments[0] - oo252*moments[1]
      + oo40*moments[2] - oo40*moments[3] - oo36*moments[5]
      - oo72*moments[6] - oo12*moments[7] - oo24*moments[8]
      + 0.25*moments[9] + 0.125*moments[12] + 0.125*moments[13];
    // Momentum to velocity
    mac[1] *= rhoInv;
    mac[2] *= rhoInv;
    mac[3] *= rhoInv;
  }
private:
  Real omega;
  static constexpr PetscScalar s1 = Constants::s1;
  static constexpr PetscScalar s2 = Constants::s2;
  static constexpr PetscScalar s4 = Constants::s4;
  static constexpr PetscScalar s10 = Constants::s10;
  static constexpr PetscScalar s16 = Constants::s16;
  static constexpr PetscScalar we = Constants::we;
  static constexpr PetscScalar wej = Constants::wej;
  static constexpr PetscScalar wxx = Constants::wxx;
  // Fractions
  static constexpr PetscScalar fo399 = 5./399.;
  static constexpr PetscScalar eo2394 = 11./2394.;
  static constexpr PetscScalar fo1197 = 4./1197.;
  static constexpr PetscScalar oo21 = 1./21.;
  static constexpr PetscScalar oo63 = 1./63.;
  static constexpr PetscScalar oo36 = 1./36.;
  static constexpr PetscScalar oo12 = 1./12.;
  static constexpr PetscScalar oo18 = 1./18.;
  static constexpr PetscScalar oo40 = 1./40.;
  static constexpr PetscScalar oo252 = 1./252.;
  static constexpr PetscScalar oo72 = 1./72.;
  static constexpr PetscScalar oo24 = 1./24.;
};

#endif
