#ifndef THERMALMRT2D
#define THERMALMRT2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/codiheader.hh"
#include "dynamics/thermalMacros2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "adjointDynamics/reverseThermalMRT2d.hh"

template <class Lattice, class Real = PetscScalar>
class ThermalMRT2d {};

template <class Real>
class ThermalMRT2d<D2Q5_TMRT,Real> {

public:
  using LatticeType = D2Q5_TMRT;
  using RealType = Real;
  using Macros = StandardThermalMacros2d<D2Q5_TMRT,Real>;
  using Equilibrium = ThermalMRTFeq2d<D2Q5_TMRT,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  ThermalMRT2d(ThermalFlowParameters par){
    Real nu = par.incPar.velocityChar * par.incPar.lengthChar /
      par.incPar.ReynoldsNumber;
    diffusivity = nu / par.PrandtlNumber;
    thermalOmega = 1./(D2Q5_TMRT::csSqInv*diffusivity + 0.5);
  }
  ThermalMRT2d(Real _om, Real _diff) : thermalOmega(_om), diffusivity(_diff){}
  ThermalMRT2d() : thermalOmega(1.), diffusivity(1.){}
  ThermalMRT2d<D2Q5_TMRT,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return ThermalMRT2d<D2Q5_TMRT,CodiReverseType<Real>>
      (this->thermalOmega,this->diffusivity);
  }
  ReverseThermalMRT2d<D2Q5_TMRT,Real> getSourceAdjoint() const
  {
    return ReverseThermalMRT2d<D2Q5_TMRT,Real>(this->thermalOmega);
  }

  Real getDiffusivity() const
  {
    return diffusivity;
  }
  Real getOmega() const
  {
    return thermalOmega;
  }
  void operator()(const Real* fdistIn, Real* fdistOut,
                  PETSC_RESTRICT Real* mac)
  {
    Macros::compute(fdistIn,mac[3]);
    Real moments[D2Q5_TMRT::numDOF - 1]; // First moment is temperature
    // Compute moments
    moments[0] = -fdistIn[2] + fdistIn[4];
    moments[1] = -fdistIn[1] + fdistIn[3];
    moments[2] = -4.*fdistIn[0] + fdistIn[1] + fdistIn[2] + fdistIn[3] + fdistIn[4];
    moments[3] = -fdistIn[1] + fdistIn[2] - fdistIn[3] + fdistIn[4];
    // Collision in moment space
    moments[0] = thermalOmega*(moments[0] - mac[1]*mac[3]);
    moments[1] = thermalOmega*(moments[1] - mac[2]*mac[3]);
    moments[2] = s3*(moments[2] + 2.*mac[3]);
    moments[3] *= s4;
    // Map back to distribution space
    fdistOut[0] = fdistIn[0] + 0.2*moments[2];
    fdistOut[1] = fdistIn[1] + 0.5*moments[1] - 0.05*moments[2]
      + 0.25*moments[3];
    fdistOut[2] = fdistIn[2] + 0.5*moments[0] - 0.05*moments[2]
      - 0.25*moments[3];
    fdistOut[3] = fdistIn[3] - 0.5*moments[1] - 0.05*moments[2]
      + 0.25*moments[3];
    fdistOut[4] = fdistIn[4] - 0.5*moments[0] - 0.05*moments[2]
      - 0.25*moments[3];
  }
  void operator()(const Real* fdistIn, Real* fdistOut,
                  PETSC_RESTRICT Real* mac, const Real& inputOmega)
  {
    Macros::compute(fdistIn,mac[3]);
    Real moments[D2Q5_TMRT::numDOF - 1]; // First moment is temperature
    // Compute moments
    moments[0] = -fdistIn[2] + fdistIn[4];
    moments[1] = -fdistIn[1] + fdistIn[3];
    moments[2] = -4.*fdistIn[0] + fdistIn[1] + fdistIn[2] + fdistIn[3] + fdistIn[4];
    moments[3] = -fdistIn[1] + fdistIn[2] - fdistIn[3] + fdistIn[4];
    // Collision in moment space
    moments[0] = inputOmega*(moments[0] - mac[1]*mac[3]);
    moments[1] = inputOmega*(moments[1] - mac[2]*mac[3]);
    moments[2] = s3*(moments[2] + 2.*mac[3]);
    moments[3] *= s4;
    // Map back to distribution space
    fdistOut[0] = fdistIn[0] + 0.2*moments[2];
    fdistOut[1] = fdistIn[1] + 0.5*moments[1] - 0.05*moments[2]
      + 0.25*moments[3];
    fdistOut[2] = fdistIn[2] + 0.5*moments[0] - 0.05*moments[2]
      - 0.25*moments[3];
    fdistOut[3] = fdistIn[3] - 0.5*moments[1] - 0.05*moments[2]
      + 0.25*moments[3];
    fdistOut[4] = fdistIn[4] - 0.5*moments[0] - 0.05*moments[2]
      - 0.25*moments[3];
  }
private:
  Real thermalOmega;
  Real diffusivity;
  static constexpr PetscScalar s3 = 1.5;
  static constexpr PetscScalar s4 = 1.5;
};

#endif
