#ifndef REVERSETHERMALBGK2D
#define REVERSETHERMALBGK2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "dynamics/thermalMacros2d.hh"
#include "dynamics/thermalEquilibriumDistributions2d.hh"
#include "adjointDynamics/reverseThermalMacros2d.hh"
#include "adjointDynamics/reverseThermalEquilibriumDistributions2d.hh"

template <class Lattice, class Real = PetscScalar>
class ReverseThermalBGK2d {

public:
  using RevMacros = ReverseThermalMacros2d<Lattice,Real>;
  using RevEquilibrium = ReverseThermalFeq2d<Lattice,Real>;
  using Macros = StandardThermalMacros2d<Lattice,Real>;
  using Equilibrium = ThermalFeq2d<Lattice,Real>;
  using RealType = Real;

  ReverseThermalBGK2d(ThermalFlowParameters par){
    Real nu = par.incPar.velocityChar * par.incPar.lengthChar /
      par.incPar.ReynoldsNumber;
    Real diffusivity = nu / par.PrandtlNumber;
    thermalOmega = 1./(Lattice::csSqInv*diffusivity + 0.5);
  }
  ReverseThermalBGK2d(Real _om) : thermalOmega(_om){}

  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         Real* mac, Real* fdistInRev,
                         Real* fdistOutRev, Real* macRev) const
  {
    Macros::compute(fdistIn,mac[3]);
    Real eqRev;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      eqRev = thermalOmega*fdistOutRev[ii];
      RevEquilibrium::eq(ii,mac[1],mac[2],mac[3],macRev[1],macRev[2],
                         macRev[3],eqRev);
      fdistInRev[ii] += (1. - thermalOmega)*fdistOutRev[ii];
    }
    RevMacros::compute(fdistInRev,macRev[3]);
  }
  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         Real* mac, Real* fdistInRev,
                         Real* fdistOutRev, Real* macRev,
                         const Real& inputOmega, Real* inputOmegaRev) const
  {
    Macros::compute(fdistIn,mac[3]);
    Real eqRev;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      *inputOmegaRev += Equilibrium::eq(ii,mac[1],mac[2],mac[3])*fdistOutRev[ii];
      eqRev = inputOmega*fdistOutRev[ii];
      RevEquilibrium::eq(ii,mac[1],mac[2],mac[3],macRev[1],macRev[2],
                         macRev[3],eqRev);
      fdistInRev[ii] += (1. - inputOmega)*fdistOutRev[ii];
      *inputOmegaRev -= fdistIn[ii]*fdistOutRev[ii];
    }
    RevMacros::compute(fdistInRev,macRev[3]);
  }
private:
  Real thermalOmega;
};

#endif
