#ifndef REVERSEMRT3D
#define REVERSEMRT3D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "dynamics/MRT3d.hh"
#include "latticeBoltzmann/lattices3d.hh"
#include "adjointDynamics/reverseIsothermalMacros3d.hh"

template <class Lattice, class Real = PetscScalar,
          class Constants = MRT3d_Compile_constants>
class ReverseIncompressibleMRT3d {};

template <class Real, class Constants>
class ReverseIncompressibleMRT3d<D3Q19,Real,Constants> {

public:
  using RevMacros = ReverseIncompressibleMacros3d<D3Q19,Real>;
  using Macros = IncompressibleMacros3d<D3Q19,Real>;
  using RealType = Real;

  ReverseIncompressibleMRT3d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  ReverseIncompressibleMRT3d(Real _om) : omega(_om){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Real tempb;
    Real tempb0;
    Real tempb1;
    Real tempb2;
    Real tempb3;
    Macros::compute(fdistIn,mac[0],mac[1],mac[2],mac[3]);
    Real momentsRev[D3Q19::numDOF - 4];
    Real jSqrb;
    Real tempEqb;
    for (size_t ii = 0; ii < D3Q19::numDOF-4; ++ii){
      momentsRev[ii] = 0.;
    }
    // Reverse map back to distributions space
    momentsRev[0] = fo399*fdistOutRev[0] + eo2394*fdistOutRev[1] + eo2394*fdistOutRev[2] +
      eo2394*fdistOutRev[3] - fo1197*fdistOutRev[4] - fo1197*fdistOutRev[5] -
      fo1197*fdistOutRev[6] - fo1197*fdistOutRev[7] - fo1197*fdistOutRev[8] -
      fo1197*fdistOutRev[9] + eo2394*fdistOutRev[10] + eo2394*fdistOutRev[11] +
      eo2394*fdistOutRev[12] - fo1197*fdistOutRev[13] - fo1197*fdistOutRev[14] -
      fo1197*fdistOutRev[15] - fo1197*fdistOutRev[16] - fo1197*fdistOutRev[17] -
      fo1197*fdistOutRev[18];
    momentsRev[1] = -oo21*fdistOutRev[0] + oo63*fdistOutRev[1] + oo63*fdistOutRev[2] +
      oo63*fdistOutRev[3] - oo252*fdistOutRev[4] - oo252*fdistOutRev[5] -
      oo252*fdistOutRev[6] - oo252*fdistOutRev[7] - oo252*fdistOutRev[8] -
      oo252*fdistOutRev[9] + oo63*fdistOutRev[10] + oo63*fdistOutRev[11] +
      oo63*fdistOutRev[12] - oo252*fdistOutRev[13] - oo252*fdistOutRev[14] -
      oo252*fdistOutRev[15] - oo252*fdistOutRev[16] - oo252*fdistOutRev[17] -
      oo252*fdistOutRev[18];
    momentsRev[2] = -0.1*fdistOutRev[3] + oo40*fdistOutRev[6] - oo40*fdistOutRev[7] +
      oo40*fdistOutRev[8] - oo40*fdistOutRev[9] + 0.1*fdistOutRev[12] -
      oo40*fdistOutRev[15] + oo40*fdistOutRev[16] - oo40*fdistOutRev[17] +
      oo40*fdistOutRev[18];
    momentsRev[3] = -0.1*fdistOutRev[2] + oo40*fdistOutRev[4] - oo40*fdistOutRev[5] +
      oo40*fdistOutRev[8] + oo40*fdistOutRev[9] + 0.1*fdistOutRev[11] -
      oo40*fdistOutRev[13] + oo40*fdistOutRev[14] - oo40*fdistOutRev[17] -
      oo40*fdistOutRev[18];
    momentsRev[4] = -0.1*fdistOutRev[1] + oo40*fdistOutRev[4] + oo40*fdistOutRev[5]
      + oo40*fdistOutRev[6] + oo40*fdistOutRev[7] + 0.1*fdistOutRev[10] -
      oo40*fdistOutRev[13] - oo40*fdistOutRev[14] - oo40*fdistOutRev[15] -
      oo40*fdistOutRev[16];
    momentsRev[5] = oo36*fdistOutRev[1] + oo36*fdistOutRev[2] - oo18*fdistOutRev[3] +
      oo18*fdistOutRev[4] + oo18*fdistOutRev[5] - oo36*fdistOutRev[6] -
      oo36*fdistOutRev[7] - oo36*fdistOutRev[8] - oo36*fdistOutRev[9] +
      oo36*fdistOutRev[10] + oo36*fdistOutRev[11] - oo18*fdistOutRev[12] +
      oo18*fdistOutRev[13] + oo18*fdistOutRev[14] - oo36*fdistOutRev[15] -
      oo36*fdistOutRev[16] - oo36*fdistOutRev[17] - oo36*fdistOutRev[18];
    momentsRev[6] = -oo36*fdistOutRev[1] - oo36*fdistOutRev[2] +
      oo18*fdistOutRev[3] + oo36*fdistOutRev[4] + oo36*fdistOutRev[5] -
      oo72*fdistOutRev[6] - oo72*fdistOutRev[7] - oo72*fdistOutRev[8] -
      oo72*fdistOutRev[9] - oo36*fdistOutRev[10] - oo36*fdistOutRev[11] +
      oo18*fdistOutRev[12] + oo36*fdistOutRev[13] + oo36*fdistOutRev[14] -
      oo72*fdistOutRev[15] - oo72*fdistOutRev[16] - oo72*fdistOutRev[17] -
      oo72*fdistOutRev[18];
    momentsRev[7] = oo12*fdistOutRev[1] - oo12*fdistOutRev[2] + oo12*fdistOutRev[6] +
      oo12*fdistOutRev[7] - oo12*fdistOutRev[8] - oo12*fdistOutRev[9] +
      oo12*fdistOutRev[10] - oo12*fdistOutRev[11] + oo12*fdistOutRev[15] +
      oo12*fdistOutRev[16] - oo12*fdistOutRev[17] - oo12*fdistOutRev[18];
    fdistInRev[18] = fdistInRev[18] + fdistOutRev[18];
    momentsRev[9] = momentsRev[9] + 0.25*fdistOutRev[18];
    momentsRev[8] = momentsRev[8] - oo24*fdistOutRev[18];
    momentsRev[12] = momentsRev[12] + 0.125*fdistOutRev[18];
    momentsRev[13] = momentsRev[13] + 0.125*fdistOutRev[18];
    fdistInRev[17] = fdistInRev[17] + fdistOutRev[17];
    momentsRev[13] = momentsRev[13] + 0.125*fdistOutRev[17];
    momentsRev[8] = momentsRev[8] - oo24*fdistOutRev[17];
    momentsRev[12] = momentsRev[12] - 0.125*fdistOutRev[17];
    momentsRev[9] = momentsRev[9] - 0.25*fdistOutRev[17];
    fdistInRev[16] = fdistInRev[16] + fdistOutRev[16];
    momentsRev[8] = momentsRev[8] + oo24*fdistOutRev[16];
    momentsRev[11] = momentsRev[11] + 0.25*fdistOutRev[16];
    momentsRev[12] = momentsRev[12] - 0.125*fdistOutRev[16];
    momentsRev[14] = momentsRev[14] - 0.125*fdistOutRev[16];
    fdistInRev[15] = fdistInRev[15] + fdistOutRev[15];
    momentsRev[8] = momentsRev[8] + oo24*fdistOutRev[15];
    momentsRev[12] = momentsRev[12] + 0.125*fdistOutRev[15];
    momentsRev[11] = momentsRev[11] - 0.25*fdistOutRev[15];
    momentsRev[14] = momentsRev[14] - 0.125*fdistOutRev[15];
    fdistInRev[14] = fdistInRev[14] + fdistOutRev[14];
    momentsRev[10] = momentsRev[10] + 0.25*fdistOutRev[14];
    momentsRev[13] = momentsRev[13] + 0.125*fdistOutRev[14];
    momentsRev[14] = momentsRev[14] + 0.125*fdistOutRev[14];
    fdistInRev[13] = fdistInRev[13] + fdistOutRev[13];
    momentsRev[10] = momentsRev[10] - 0.25*fdistOutRev[13];
    momentsRev[14] = momentsRev[14] + 0.125*fdistOutRev[13];
    momentsRev[13] = momentsRev[13] - 0.125*fdistOutRev[13];
    fdistInRev[12] = fdistInRev[12] + fdistOutRev[12];
    fdistInRev[11] = fdistInRev[11] + fdistOutRev[11];
    momentsRev[8] = momentsRev[8] + oo12*fdistOutRev[11];
    fdistInRev[10] = fdistInRev[10] + fdistOutRev[10];
    momentsRev[8] = momentsRev[8] - oo12*fdistOutRev[10];
    fdistInRev[9] = fdistInRev[9] + fdistOutRev[9];
    momentsRev[9] = momentsRev[9] + 0.25*fdistOutRev[9];
    momentsRev[8] = momentsRev[8] - oo24*fdistOutRev[9];
    momentsRev[12] = momentsRev[12] - 0.125*fdistOutRev[9];
    momentsRev[13] = momentsRev[13] - 0.125*fdistOutRev[9];
    fdistInRev[8] = fdistInRev[8] + fdistOutRev[8];
    momentsRev[12] = momentsRev[12] + 0.125*fdistOutRev[8];
    momentsRev[8] = momentsRev[8] - oo24*fdistOutRev[8];
    momentsRev[9] = momentsRev[9] - 0.25*fdistOutRev[8];
    momentsRev[13] = momentsRev[13] - 0.125*fdistOutRev[8];
    fdistInRev[7] = fdistInRev[7] + fdistOutRev[7];
    momentsRev[8] = momentsRev[8] + oo24*fdistOutRev[7];
    momentsRev[11] = momentsRev[11] + 0.25*fdistOutRev[7];
    momentsRev[12] = momentsRev[12] + 0.125*fdistOutRev[7];
    momentsRev[14] = momentsRev[14] + 0.125*fdistOutRev[7];
    fdistInRev[6] = fdistInRev[6] + fdistOutRev[6];
    momentsRev[8] = momentsRev[8] + oo24*fdistOutRev[6];
    momentsRev[14] = momentsRev[14] + 0.125*fdistOutRev[6];
    momentsRev[12] = momentsRev[12] - 0.125*fdistOutRev[6];
    momentsRev[11] = momentsRev[11] - 0.25*fdistOutRev[6];
    fdistInRev[5] = fdistInRev[5] + fdistOutRev[5];
    momentsRev[10] = momentsRev[10] + 0.25*fdistOutRev[5];
    momentsRev[13] = momentsRev[13] - 0.125*fdistOutRev[5];
    momentsRev[14] = momentsRev[14] - 0.125*fdistOutRev[5];
    fdistInRev[4] = fdistInRev[4] + fdistOutRev[4];
    momentsRev[10] = momentsRev[10] - 0.25*fdistOutRev[4];
    momentsRev[13] = momentsRev[13] + 0.125*fdistOutRev[4];
    momentsRev[14] = momentsRev[14] - 0.125*fdistOutRev[4];
    fdistInRev[3] = fdistInRev[3] + fdistOutRev[3];
    fdistInRev[2] = fdistInRev[2] + fdistOutRev[2];
    momentsRev[8] = momentsRev[8] + oo12*fdistOutRev[2];
    fdistInRev[1] = fdistInRev[1] + fdistOutRev[1];
    momentsRev[8] = momentsRev[8] - oo12*fdistOutRev[1];
    fdistInRev[0] = fdistInRev[0] + fdistOutRev[0];
    momentsRev[14] = s16*momentsRev[14];
    momentsRev[13] = s16*momentsRev[13];
    momentsRev[12] = s16*momentsRev[12];
    tempb = omega*momentsRev[11];
    macRev[1] = macRev[1] - mac[3]*tempb;
    macRev[3] = macRev[3] - mac[1]*tempb;
    momentsRev[11] = tempb;
    tempb0 = omega*momentsRev[10];
    macRev[2] = macRev[2] - mac[3]*tempb0;
    macRev[3] = macRev[3] - mac[2]*tempb0;
    momentsRev[10] = tempb0;
    tempb1 = omega*momentsRev[9];
    macRev[1] = macRev[1] - mac[2]*tempb1;
    macRev[2] = macRev[2] - mac[1]*tempb1;
    momentsRev[9] = tempb1;
    tempEqb = -(s10*wxx*momentsRev[8]);
    momentsRev[8] = s10*momentsRev[8];
    tempEqb = tempEqb - omega*momentsRev[7];
    momentsRev[7] = omega*momentsRev[7];
    macRev[2] = macRev[2] + 2*mac[2]*tempEqb;
    macRev[3] = macRev[3] - 2*mac[3]*tempEqb;
    tempEqb = -(s10*wxx*momentsRev[6]);
    momentsRev[6] = s10*momentsRev[6];
    tempEqb = tempEqb - omega*momentsRev[5];
    momentsRev[5] = omega*momentsRev[5];
    macRev[1] = macRev[1] + 2.*2*mac[1]*tempEqb;
    macRev[2] = macRev[2] - 2*mac[2]*tempEqb;
    macRev[3] = macRev[3] - 2*mac[3]*tempEqb;
    macRev[3] = macRev[3] + 2.*s4*momentsRev[4]/3.;
    momentsRev[4] = s4*momentsRev[4];
    macRev[2] = macRev[2] + 2.*s4*momentsRev[3]/3.;
    momentsRev[3] = s4*momentsRev[3];
    macRev[1] = macRev[1] + 2.*s4*momentsRev[2]/3.;
    momentsRev[2] = s4*momentsRev[2];
    tempb2 = s2*momentsRev[1];
    macRev[0] = macRev[0] - we*tempb2;
    momentsRev[1] = tempb2;
    tempb3 = s1*momentsRev[0];
    jSqrb = -(19.*tempb3) - wej*tempb2;
    macRev[0] = macRev[0] + 11.*tempb3;
    momentsRev[0] = tempb3;
    macRev[1] = macRev[1] + 2*mac[1]*jSqrb;
    macRev[2] = macRev[2] + 2*mac[2]*jSqrb;
    macRev[3] = macRev[3] + 2*mac[3]*jSqrb;
    fdistInRev[4] = fdistInRev[4] + momentsRev[14];
    fdistInRev[5] = fdistInRev[5] + momentsRev[14];
    fdistInRev[6] = fdistInRev[6] - momentsRev[14];
    fdistInRev[7] = fdistInRev[7] - momentsRev[14];
    fdistInRev[15] = fdistInRev[15] + momentsRev[14];
    fdistInRev[14] = fdistInRev[14] - momentsRev[14];
    fdistInRev[13] = fdistInRev[13] - momentsRev[14];
    fdistInRev[16] = fdistInRev[16] + momentsRev[14];
    fdistInRev[5] = fdistInRev[5] + momentsRev[13];
    fdistInRev[4] = fdistInRev[4] - momentsRev[13];
    fdistInRev[8] = fdistInRev[8] + momentsRev[13];
    fdistInRev[9] = fdistInRev[9] + momentsRev[13];
    fdistInRev[13] = fdistInRev[13] + momentsRev[13];
    fdistInRev[14] = fdistInRev[14] - momentsRev[13];
    fdistInRev[17] = fdistInRev[17] - momentsRev[13];
    fdistInRev[18] = fdistInRev[18] - momentsRev[13];
    fdistInRev[6] = fdistInRev[6] + momentsRev[12];
    fdistInRev[7] = fdistInRev[7] - momentsRev[12];
    fdistInRev[9] = fdistInRev[9] + momentsRev[12];
    fdistInRev[8] = fdistInRev[8] - momentsRev[12];
    fdistInRev[16] = fdistInRev[16] + momentsRev[12];
    fdistInRev[15] = fdistInRev[15] - momentsRev[12];
    fdistInRev[17] = fdistInRev[17] + momentsRev[12];
    fdistInRev[18] = fdistInRev[18] - momentsRev[12];
    fdistInRev[6] = fdistInRev[6] + momentsRev[11];
    fdistInRev[7] = fdistInRev[7] - momentsRev[11];
    fdistInRev[15] = fdistInRev[15] + momentsRev[11];
    fdistInRev[16] = fdistInRev[16] - momentsRev[11];
    fdistInRev[4] = fdistInRev[4] + momentsRev[10];
    fdistInRev[5] = fdistInRev[5] - momentsRev[10];
    fdistInRev[13] = fdistInRev[13] + momentsRev[10];
    fdistInRev[14] = fdistInRev[14] - momentsRev[10];
    fdistInRev[8] = fdistInRev[8] + momentsRev[9];
    fdistInRev[9] = fdistInRev[9] - momentsRev[9];
    fdistInRev[17] = fdistInRev[17] + momentsRev[9];
    fdistInRev[18] = fdistInRev[18] - momentsRev[9];
    fdistInRev[1] = fdistInRev[1] + 2.*momentsRev[8];
    fdistInRev[2] = fdistInRev[2] - 2.*momentsRev[8];
    fdistInRev[6] = fdistInRev[6] - momentsRev[8];
    fdistInRev[8] = fdistInRev[8] + momentsRev[8];
    fdistInRev[7] = fdistInRev[7] - momentsRev[8];
    fdistInRev[9] = fdistInRev[9] + momentsRev[8];
    fdistInRev[10] = fdistInRev[10] + 2.*momentsRev[8];
    fdistInRev[11] = fdistInRev[11] - 2.*momentsRev[8];
    fdistInRev[15] = fdistInRev[15] - momentsRev[8];
    fdistInRev[17] = fdistInRev[17] + momentsRev[8];
    fdistInRev[16] = fdistInRev[16] - momentsRev[8];
    fdistInRev[18] = fdistInRev[18] + momentsRev[8];
    fdistInRev[2] = fdistInRev[2] + momentsRev[7];
    fdistInRev[1] = fdistInRev[1] - momentsRev[7];
    fdistInRev[6] = fdistInRev[6] - momentsRev[7];
    fdistInRev[8] = fdistInRev[8] + momentsRev[7];
    fdistInRev[7] = fdistInRev[7] - momentsRev[7];
    fdistInRev[9] = fdistInRev[9] + momentsRev[7];
    fdistInRev[11] = fdistInRev[11] + momentsRev[7];
    fdistInRev[10] = fdistInRev[10] - momentsRev[7];
    fdistInRev[15] = fdistInRev[15] - momentsRev[7];
    fdistInRev[17] = fdistInRev[17] + momentsRev[7];
    fdistInRev[16] = fdistInRev[16] - momentsRev[7];
    fdistInRev[18] = fdistInRev[18] + momentsRev[7];
    fdistInRev[1] = fdistInRev[1] + 2.*momentsRev[6];
    fdistInRev[2] = fdistInRev[2] + 2.*momentsRev[6];
    fdistInRev[3] = fdistInRev[3] - 4.*momentsRev[6];
    fdistInRev[4] = fdistInRev[4] - 2.*momentsRev[6];
    fdistInRev[5] = fdistInRev[5] - 2.*momentsRev[6];
    fdistInRev[6] = fdistInRev[6] + momentsRev[6];
    fdistInRev[7] = fdistInRev[7] + momentsRev[6];
    fdistInRev[8] = fdistInRev[8] + momentsRev[6];
    fdistInRev[9] = fdistInRev[9] + momentsRev[6];
    fdistInRev[10] = fdistInRev[10] + 2.*momentsRev[6];
    fdistInRev[11] = fdistInRev[11] + 2.*momentsRev[6];
    fdistInRev[12] = fdistInRev[12] - 4.*momentsRev[6];
    fdistInRev[13] = fdistInRev[13] - 2.*momentsRev[6];
    fdistInRev[14] = fdistInRev[14] - 2.*momentsRev[6];
    fdistInRev[15] = fdistInRev[15] + momentsRev[6];
    fdistInRev[16] = fdistInRev[16] + momentsRev[6];
    fdistInRev[17] = fdistInRev[17] + momentsRev[6];
    fdistInRev[18] = fdistInRev[18] + momentsRev[6];
    fdistInRev[3] = fdistInRev[3] + 2.*momentsRev[5];
    fdistInRev[2] = fdistInRev[2] - momentsRev[5];
    fdistInRev[1] = fdistInRev[1] - momentsRev[5];
    fdistInRev[4] = fdistInRev[4] - 2.*momentsRev[5];
    fdistInRev[5] = fdistInRev[5] - 2.*momentsRev[5];
    fdistInRev[6] = fdistInRev[6] + momentsRev[5];
    fdistInRev[7] = fdistInRev[7] + momentsRev[5];
    fdistInRev[8] = fdistInRev[8] + momentsRev[5];
    fdistInRev[9] = fdistInRev[9] + momentsRev[5];
    fdistInRev[12] = fdistInRev[12] + 2.*momentsRev[5];
    fdistInRev[11] = fdistInRev[11] - momentsRev[5];
    fdistInRev[10] = fdistInRev[10] - momentsRev[5];
    fdistInRev[13] = fdistInRev[13] - 2.*momentsRev[5];
    fdistInRev[14] = fdistInRev[14] - 2.*momentsRev[5];
    fdistInRev[15] = fdistInRev[15] + momentsRev[5];
    fdistInRev[16] = fdistInRev[16] + momentsRev[5];
    fdistInRev[17] = fdistInRev[17] + momentsRev[5];
    fdistInRev[18] = fdistInRev[18] + momentsRev[5];
    fdistInRev[1] = fdistInRev[1] + 4.*momentsRev[4];
    fdistInRev[4] = fdistInRev[4] - momentsRev[4];
    fdistInRev[5] = fdistInRev[5] - momentsRev[4];
    fdistInRev[6] = fdistInRev[6] - momentsRev[4];
    fdistInRev[7] = fdistInRev[7] - momentsRev[4];
    fdistInRev[13] = fdistInRev[13] + momentsRev[4];
    fdistInRev[10] = fdistInRev[10] - 4.*momentsRev[4];
    fdistInRev[14] = fdistInRev[14] + momentsRev[4];
    fdistInRev[15] = fdistInRev[15] + momentsRev[4];
    fdistInRev[16] = fdistInRev[16] + momentsRev[4];
    fdistInRev[2] = fdistInRev[2] + 4.*momentsRev[3];
    fdistInRev[4] = fdistInRev[4] - momentsRev[3];
    fdistInRev[5] = fdistInRev[5] + momentsRev[3];
    fdistInRev[8] = fdistInRev[8] - momentsRev[3];
    fdistInRev[9] = fdistInRev[9] - momentsRev[3];
    fdistInRev[13] = fdistInRev[13] + momentsRev[3];
    fdistInRev[11] = fdistInRev[11] - 4.*momentsRev[3];
    fdistInRev[14] = fdistInRev[14] - momentsRev[3];
    fdistInRev[17] = fdistInRev[17] + momentsRev[3];
    fdistInRev[18] = fdistInRev[18] + momentsRev[3];
    fdistInRev[3] = fdistInRev[3] + 4.*momentsRev[2];
    fdistInRev[6] = fdistInRev[6] - momentsRev[2];
    fdistInRev[7] = fdistInRev[7] + momentsRev[2];
    fdistInRev[9] = fdistInRev[9] + momentsRev[2];
    fdistInRev[8] = fdistInRev[8] - momentsRev[2];
    fdistInRev[15] = fdistInRev[15] + momentsRev[2];
    fdistInRev[12] = fdistInRev[12] - 4.*momentsRev[2];
    fdistInRev[16] = fdistInRev[16] - momentsRev[2];
    fdistInRev[17] = fdistInRev[17] + momentsRev[2];
    fdistInRev[18] = fdistInRev[18] - momentsRev[2];
    fdistInRev[0] = fdistInRev[0] + 12.*momentsRev[1];
    fdistInRev[1] = fdistInRev[1] - 4.*momentsRev[1];
    fdistInRev[2] = fdistInRev[2] - 4.*momentsRev[1];
    fdistInRev[4] = fdistInRev[4] + momentsRev[1];
    fdistInRev[3] = fdistInRev[3] - 4.*momentsRev[1];
    fdistInRev[5] = fdistInRev[5] + momentsRev[1];
    fdistInRev[6] = fdistInRev[6] + momentsRev[1];
    fdistInRev[7] = fdistInRev[7] + momentsRev[1];
    fdistInRev[8] = fdistInRev[8] + momentsRev[1];
    fdistInRev[9] = fdistInRev[9] + momentsRev[1];
    fdistInRev[13] = fdistInRev[13] + momentsRev[1];
    fdistInRev[11] = fdistInRev[11] - 4.*momentsRev[1];
    fdistInRev[12] = fdistInRev[12] - 4.*momentsRev[1];
    fdistInRev[10] = fdistInRev[10] - 4.*momentsRev[1];
    fdistInRev[14] = fdistInRev[14] + momentsRev[1];
    fdistInRev[15] = fdistInRev[15] + momentsRev[1];
    fdistInRev[16] = fdistInRev[16] + momentsRev[1];
    fdistInRev[17] = fdistInRev[17] + momentsRev[1];
    fdistInRev[18] = fdistInRev[18] + momentsRev[1];
    fdistInRev[4] = fdistInRev[4] + 8.*momentsRev[0];
    fdistInRev[1] = fdistInRev[1] - 11.*momentsRev[0];
    fdistInRev[2] = fdistInRev[2] - 11.*momentsRev[0];
    fdistInRev[0] = fdistInRev[0] - 30.*momentsRev[0];
    fdistInRev[3] = fdistInRev[3] - 11.*momentsRev[0];
    fdistInRev[5] = fdistInRev[5] + 8.*momentsRev[0];
    fdistInRev[6] = fdistInRev[6] + 8.*momentsRev[0];
    fdistInRev[7] = fdistInRev[7] + 8.*momentsRev[0];
    fdistInRev[8] = fdistInRev[8] + 8.*momentsRev[0];
    fdistInRev[9] = fdistInRev[9] + 8.*momentsRev[0];
    fdistInRev[13] = fdistInRev[13] + 8.*momentsRev[0];
    fdistInRev[11] = fdistInRev[11] - 11.*momentsRev[0];
    fdistInRev[12] = fdistInRev[12] - 11.*momentsRev[0];
    fdistInRev[10] = fdistInRev[10] - 11.*momentsRev[0];
    fdistInRev[14] = fdistInRev[14] + 8.*momentsRev[0];
    fdistInRev[15] = fdistInRev[15] + 8.*momentsRev[0];
    fdistInRev[16] = fdistInRev[16] + 8.*momentsRev[0];
    fdistInRev[17] = fdistInRev[17] + 8.*momentsRev[0];
    fdistInRev[18] = fdistInRev[18] + 8.*momentsRev[0];
    RevMacros::compute(fdistInRev,macRev[0],macRev[1],macRev[2],macRev[3]);
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
