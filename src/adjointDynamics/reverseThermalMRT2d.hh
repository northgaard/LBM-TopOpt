#ifndef REVERSETHERMALMRT2D
#define REVERSETHERMALMRT2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "dynamics/thermalMacros2d.hh"
#include "adjointDynamics/reverseThermalMacros2d.hh"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
class ReverseThermalMRT2d {};

template <class Real>
class ReverseThermalMRT2d<D2Q5_TMRT,Real> {

public:
  using RevMacros = ReverseThermalMacros2d<D2Q5_TMRT,Real>;
  using Macros = StandardThermalMacros2d<D2Q5_TMRT,Real>;
  using RealType = Real;

  ReverseThermalMRT2d(ThermalFlowParameters par){
    Real nu = par.incPar.velocityChar * par.incPar.lengthChar /
      par.incPar.ReynoldsNumber;
    Real diffusivity = nu / par.PrandtlNumber;
    thermalOmega = 1./(D2Q5_TMRT::csSqInv*diffusivity + 0.5);
  }
  ReverseThermalMRT2d(Real _om) : thermalOmega(_om){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Macros::compute(fdistIn,mac[3]);
    Real momentsRev[D2Q5_TMRT::numDOF - 1];
    // Reverse map back to distribution space
    momentsRev[0] = 0.5*fdistOutRev[2] - 0.5*fdistOutRev[4];
    momentsRev[1] = 0.5*fdistOutRev[1] - 0.5*fdistOutRev[3];
    momentsRev[2] = 0.2*fdistOutRev[0] - 0.05*fdistOutRev[1] - 0.05*fdistOutRev[2]
      - 0.05*fdistOutRev[3] - 0.05*fdistOutRev[4];
    momentsRev[3] = 0.25*(fdistOutRev[1] - fdistOutRev[2] + fdistOutRev[3]
                          - fdistOutRev[4]);
    // Reverse collision in moment space
    momentsRev[3] *= s4;
    momentsRev[2] *= s3;
    momentsRev[1] *= thermalOmega;
    momentsRev[0] *= thermalOmega;
    macRev[1] -= mac[3]*momentsRev[0];
    macRev[2] -= mac[3]*momentsRev[1];
    macRev[3] -= mac[1]*momentsRev[0] + mac[2]*momentsRev[1] - 2.*momentsRev[2];
    // Reverse moment transform
    fdistInRev[0] += fdistOutRev[0] - 4.*momentsRev[2];
    fdistInRev[1] += fdistOutRev[1] - momentsRev[1] + momentsRev[2] - momentsRev[3];
    fdistInRev[2] += fdistOutRev[2] - momentsRev[0] + momentsRev[2] + momentsRev[3];
    fdistInRev[3] += fdistOutRev[3] + momentsRev[1] + momentsRev[2] - momentsRev[3];
    fdistInRev[4] += fdistOutRev[4] + momentsRev[0] + momentsRev[2] + momentsRev[3];
    RevMacros::compute(fdistInRev,macRev[3]);
  }
  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev,
                  const Real& inputOmega, Real* inputOmegaRev) const
  {
    Macros::compute(fdistIn,mac[3]);
    Real momentsRev[D2Q5_TMRT::numDOF - 1];
    // Reverse map back to distribution space
    momentsRev[0] = 0.5*fdistOutRev[2] - 0.5*fdistOutRev[4];
    momentsRev[1] = 0.5*fdistOutRev[1] - 0.5*fdistOutRev[3];
    momentsRev[2] = 0.2*fdistOutRev[0] - 0.05*fdistOutRev[1] - 0.05*fdistOutRev[2]
      - 0.05*fdistOutRev[3] - 0.05*fdistOutRev[4];
    momentsRev[3] = 0.25*(fdistOutRev[1] - fdistOutRev[2] + fdistOutRev[3]
                          - fdistOutRev[4]);
    *inputOmegaRev += (fdistIn[4] - fdistIn[2] - mac[2]*mac[3])*momentsRev[0]
      + (fdistIn[3] - fdistIn[1] - mac[1]*mac[3])*momentsRev[1];
    // Reverse collision in moment space
    momentsRev[3] *= s4;
    momentsRev[2] *= s3;
    momentsRev[1] *= inputOmega;
    momentsRev[0] *= inputOmega;
    macRev[1] -= mac[3]*momentsRev[0];
    macRev[2] -= mac[3]*momentsRev[1];
    macRev[3] -= mac[1]*momentsRev[0] + mac[2]*momentsRev[1] - 2.*momentsRev[2];
    // Reverse moment transform
    fdistInRev[0] += fdistOutRev[0] - 4.*momentsRev[2];
    fdistInRev[1] += fdistOutRev[1] - momentsRev[1] + momentsRev[2] - momentsRev[3];
    fdistInRev[2] += fdistOutRev[2] - momentsRev[0] + momentsRev[2] + momentsRev[3];
    fdistInRev[3] += fdistOutRev[3] + momentsRev[1] + momentsRev[2] - momentsRev[3];
    fdistInRev[4] += fdistOutRev[4] + momentsRev[0] + momentsRev[2] + momentsRev[3];
    RevMacros::compute(fdistInRev,macRev[3]);
  }
private:
  Real thermalOmega;
  static constexpr PetscScalar s3 = 1.5;
  static constexpr PetscScalar s4 = 1.5;
};

#endif
