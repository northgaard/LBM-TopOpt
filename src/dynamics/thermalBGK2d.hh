#ifndef THERMALBGK2D
#define THERMALBGK2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/codiheader.hh"
#include "dynamics/thermalMacros2d.hh"
#include "dynamics/thermalEquilibriumDistributions2d.hh"
#include "adjointDynamics/reverseThermalBGK2d.hh"

template <class Lattice, class Real = PetscScalar>
class ThermalBGK2d {

public:
  using LatticeType = Lattice;
  using RealType = Real;
  using Equilibrium = ThermalFeq2d<Lattice,Real>;
  using Macros = StandardThermalMacros2d<Lattice,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  ThermalBGK2d(ThermalFlowParameters par){
    Real nu = par.incPar.velocityChar * par.incPar.lengthChar /
      par.incPar.ReynoldsNumber;
    diffusivity = nu / par.PrandtlNumber;
    thermalOmega = 1./(Lattice::csSqInv*diffusivity + 0.5);
  }
  ThermalBGK2d(Real _om, Real _diff) : thermalOmega(_om), diffusivity(_diff){}
  ThermalBGK2d() : thermalOmega(1.), diffusivity(1.){}
  ThermalBGK2d<Lattice,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return ThermalBGK2d<Lattice,CodiReverseType<Real>>
      (this->thermalOmega,this->diffusivity);
  }
  ReverseThermalBGK2d<Lattice,Real> getSourceAdjoint() const
  {
    return ReverseThermalBGK2d<Lattice,Real>(this->thermalOmega);
  }

  inline Real getDiffusivity() const
  {
    return diffusivity;
  }
  inline Real getOmega() const
  {
    return thermalOmega;
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
			 Real* mac) const
  {
    Macros::compute(fdistIn,mac[3]);
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistOut[ii] = fdistIn[ii] * (1. - thermalOmega);
      fdistOut[ii] += thermalOmega * Equilibrium::eq(ii,mac[1],mac[2],mac[3]);
    }
  }
  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         Real* mac, const Real& inputOmega) const
  {
    Macros::compute(fdistIn,mac[3]);
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistOut[ii] = fdistIn[ii] * (1. - inputOmega);
      fdistOut[ii] += inputOmega * Equilibrium::eq(ii,mac[1],mac[2],mac[3]);
    }
  }
private:
  Real thermalOmega;
  Real diffusivity;
};

#endif
