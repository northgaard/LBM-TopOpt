#ifndef REVERSEINCOMPRESSIBLEBGK2D
#define REVERSEINCOMPRESSIBLEBGK2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "adjointDynamics/reverseEquilibriumDistributions2d.hh"
#include "adjointDynamics/reverseIsothermalMacros2d.hh"
#include "external/tapenadeHeader.hh"

template <class Lattice, class Real = PetscScalar>
class ReverseIncompressibleBGK2d {

public:
  using RevEquilibrium = ReverseIncompressibleFeq2d<Lattice,Real>;
  using RevMacros = ReverseIncompressibleMacros2d<Lattice,Real>;
  using Macros = IncompressibleMacros2d<Lattice,Real>;
  using RealType = Real;

  ReverseIncompressibleBGK2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  ReverseIncompressibleBGK2d(Real _om) : omega(_om){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqr;
    Real uSqrRev = 0.;
    Real eqRev;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      eqRev = omega*fdistOutRev[ii];
      RevEquilibrium::eq(ii,mac[0],mac[1],mac[2],uSqr,
                         macRev[0],macRev[1],macRev[2],uSqrRev,eqRev);
      fdistInRev[ii] += (1. - omega)*fdistOutRev[ii];
    }
    macRev[1] += 2.*mac[1]*uSqrRev;
    macRev[2] += 2.*mac[2]*uSqrRev;
    RevMacros::compute(fdistInRev,macRev[0],macRev[1],macRev[2]);
  }
private:
  Real omega;
};

template <class Real>
class ReverseIncompressibleBGK2d<D2Q9,Real> {

public:
  using RevEquilibrium = ReverseIncompressibleFeq2d<D2Q9,Real>;
  using RevMacros = ReverseIncompressibleMacros2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  using RealType = Real;

  ReverseIncompressibleBGK2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  ReverseIncompressibleBGK2d(Real _om) : omega(_om){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    const Real rel = 1. - omega;
    Real uSqrRev, cFeqRev;
    Real temp;

    fdistInRev[0] += rel*fdistOutRev[0];
    fdistInRev[1] += rel*fdistOutRev[1];
    fdistInRev[2] += rel*fdistOutRev[2];
    fdistInRev[3] += rel*fdistOutRev[3];
    fdistInRev[4] += rel*fdistOutRev[4];
    fdistInRev[5] += rel*fdistOutRev[5];
    fdistInRev[6] += rel*fdistOutRev[6];
    fdistInRev[7] += rel*fdistOutRev[7];
    fdistInRev[8] += rel*fdistOutRev[8];

    temp = omega*diagWeight*fdistOutRev[7];
    cFeqRev = temp;
    uSqrRev = temp;
    macRev[1] += temp*(9.*mac[2] + 3.);
    macRev[2] += 3.*temp*(3.*mac[1] + 1.);
    temp = omega*diagWeight*fdistOutRev[5];
    cFeqRev += temp;
    uSqrRev += temp;
    macRev[1] -= temp*(9.*mac[2] + 3.);
    macRev[2] += temp*(3. - 9.*mac[1]);
    temp = omega*diagWeight*fdistOutRev[3];
    cFeqRev += temp;
    uSqrRev += temp;
    macRev[1] += temp*(9.*mac[2] - 3.);
    macRev[2] += temp*(9.*mac[1] - 3.);
    temp = omega*diagWeight*fdistOutRev[1];
    cFeqRev += temp;
    uSqrRev += temp;
    uSqrRev *= 4.5;
    macRev[1] += temp*(3 - 9.*mac[2]);
    macRev[2] -= temp*(9.*mac[1] + 3.);
    temp = omega*stWeight*fdistOutRev[8];
    cFeqRev += temp;
    macRev[1] += (9.*mac[1]+3.)*temp;
    temp = omega*stWeight*fdistOutRev[6];
    cFeqRev += temp;
    macRev[2] = macRev[2] + (9.*mac[2]+3.)*temp;
    temp = omega*stWeight*fdistOutRev[4];
    cFeqRev += temp;
    macRev[1] = macRev[1] + (9.*mac[1]-3.)*temp;
    temp = omega*stWeight*fdistOutRev[2];
    cFeqRev += temp + omega*restWeight*fdistOutRev[0];
    macRev[2] = macRev[2] + (9.*mac[2]-3.)*temp;
    uSqrRev -= 1.5*cFeqRev;
    macRev[0] = macRev[0] + cFeqRev;
    macRev[1] = macRev[1] + 2.*mac[1]*uSqrRev;
    macRev[2] = macRev[2] + 2.*mac[2]*uSqrRev;
    RevMacros::compute(fdistInRev,macRev[0],macRev[1],macRev[2]);
  }
private:
  Real omega;
  static constexpr PetscScalar restWeight = 4./9.;
  static constexpr PetscScalar stWeight = 1./9.;
  static constexpr PetscScalar diagWeight = 1./36.;
};

#endif
