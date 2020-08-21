#ifndef BGK3D
#define BGK3D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "dynamics/isothermalMacros3d.hh"
#include "dynamics/equilibriumDistributions3d.hh"
#include "latticeBoltzmann/lattices3d.hh"

template <class Lattice, class Real = PetscScalar>
class StandardBGK3d {

public:
  using RealType = Real;
  using Equilibrium = StandardFeq3d<Lattice,Real>;
  using Macros = StandardMacros3d<Lattice,Real>;

  StandardBGK3d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  StandardBGK3d(Real _om) : omega(_om){}
  StandardBGK3d() : omega(1.){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac) const
  {
    Macros::compute(fdistIn,mac[0],mac[1],mac[2],mac[3]);
    Real uSqr = mac[1]*mac[1] + mac[2]*mac[2] + mac[3]*mac[3];
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistOut[ii] = fdistIn[ii] * (1. - omega);
      fdistOut[ii] += omega * Equilibrium::eq(ii,mac[0],mac[1],mac[2],mac[3],uSqr);
    }
  }
private:
  Real omega;
};

#endif
