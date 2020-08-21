#ifndef INCOMPRESSIBLEBGK2D
#define INCOMPRESSIBLEBGK2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/codiheader.hh"
#include "dynamics/isothermalMacros2d.hh"
#include "dynamics/equilibriumDistributions2d.hh"
#include "adjointDynamics/reverseIncompressibleBGK2d.hh"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
class IncompressibleBGK2d {

public:
  using LatticeType = Lattice;
  using RealType = Real;
  using Equilibrium = IncompressibleFeq2d<Lattice,Real>;
  using Macros = IncompressibleMacros2d<Lattice,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  IncompressibleBGK2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  IncompressibleBGK2d(Real _om) : omega(_om){}
  IncompressibleBGK2d() : omega(1.){}
  IncompressibleBGK2d<Lattice,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return IncompressibleBGK2d<Lattice,CodiReverseType<Real>>(this->omega);
  }
  ReverseIncompressibleBGK2d<Lattice,Real> getSourceAdjoint() const
  {
    return ReverseIncompressibleBGK2d<Lattice,Real>(this->omega);
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         PETSC_RESTRICT Real* mac) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqr = mac[1]*mac[1] + mac[2]*mac[2];
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistOut[ii] = fdistIn[ii] * (1. - omega);
      fdistOut[ii] += omega * Equilibrium::eq(ii,mac[0],mac[1],mac[2],uSqr);
    }
  }
private:
  Real omega;
};

template <class Lattice, class Real = PetscScalar>
class StandardBGK2d {

public:
  using LatticeType = Lattice;
  using RealType = Real;
  using Equilibrium = StandardFeq2d<Lattice,Real>;
  using Macros = StandardMacros2d<Lattice,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  StandardBGK2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  StandardBGK2d(Real _om) : omega(_om){}
  StandardBGK2d() : omega(1.){}
  StandardBGK2d<Lattice,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return StandardBGK2d<Lattice,CodiReverseType<Real>>(this->omega);
  }
  void getSourceAdjoint() const
  {}

  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         Real* mac) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqr = mac[1]*mac[1] + mac[2]*mac[2];
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistOut[ii] = fdistIn[ii] * (1. - omega);
      fdistOut[ii] += omega * Equilibrium::eq(ii,mac[0],mac[1],mac[2],uSqr);
    }
  }
private:
  Real omega;
};

/*
  Optimization for D2Q9
*/

template <class Real>
class IncompressibleBGK2d<D2Q9,Real> {

public:
  using LatticeType = D2Q9;
  using RealType = Real;
  using Equilibrium = IncompressibleFeq2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  static constexpr PetscInt numAdditionalFields = 0;

  IncompressibleBGK2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  IncompressibleBGK2d(Real _om) : omega(_om){}
  IncompressibleBGK2d() : omega(1.){}
  IncompressibleBGK2d<D2Q9,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return IncompressibleBGK2d<D2Q9,CodiReverseType<Real>>(this->omega);
  }
  ReverseIncompressibleBGK2d<D2Q9,Real> getSourceAdjoint() const
  {
    return ReverseIncompressibleBGK2d<D2Q9,Real>(this->omega);
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
                         PETSC_RESTRICT Real* mac) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    const Real uSqr = mac[1]*mac[1] + mac[2]*mac[2];
    const Real cFeq = mac[0] - 1.5*uSqr;
    const Real rel = 1. - omega;
    // Rest collision
    fdistOut[0] = rel*fdistIn[0] + omega*restWeight*cFeq;
    // Straight collision
    fdistOut[2] = rel*fdistIn[2] +
      omega*stWeight*(cFeq - 3.*mac[2] + 4.5*mac[2]*mac[2]);
    fdistOut[4] = rel*fdistIn[4] +
      omega*stWeight*(cFeq - 3.*mac[1] + 4.5*mac[1]*mac[1]);
    fdistOut[6] = rel*fdistIn[6] +
      omega*stWeight*(cFeq + 3.*mac[2] + 4.5*mac[2]*mac[2]);
    fdistOut[8] = rel*fdistIn[8] +
      omega*stWeight*(cFeq + 3.*mac[1] + 4.5*mac[1]*mac[1]);
    // Diagonal collision
    fdistOut[1] = rel*fdistIn[1] +
      omega*diagWeight*(cFeq + 3.*(mac[1] - mac[2]) + 4.5*(uSqr - 2.*mac[1]*mac[2]));
    fdistOut[3] = rel*fdistIn[3] +
      omega*diagWeight*(cFeq - 3.*(mac[1] + mac[2]) + 4.5*(uSqr + 2.*mac[1]*mac[2]));
    fdistOut[5] = rel*fdistIn[5] +
      omega*diagWeight*(cFeq + 3.*(mac[2] - mac[1]) + 4.5*(uSqr - 2.*mac[1]*mac[2]));
    fdistOut[7] = rel*fdistIn[7] +
      omega*diagWeight*(cFeq + 3.*(mac[1] + mac[2]) + 4.5*(uSqr + 2.*mac[1]*mac[2]));
  }
private:
  Real omega;
  static constexpr PetscScalar restWeight = 4./9.;
  static constexpr PetscScalar stWeight = 1./9.;
  static constexpr PetscScalar diagWeight = 1./36.;
};

#endif
