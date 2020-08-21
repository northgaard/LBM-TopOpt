#ifndef CASCADED2D
#define CASCADED2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/codiheader.hh"
#include "dynamics/isothermalMacros2d.hh"
#include "dynamics/equilibriumDistributions2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "adjointDynamics/reverseCascaded2d.hh"

template <class Lattice, class Real = PetscScalar>
class IncompressibleCascaded2d {};

template <class Real>
class IncompressibleCascaded2d<D2Q9,Real> {

public:
  using LatticeType = D2Q9;
  using RealType = Real;
  using Equilibrium = IncompressibleFeq2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  static constexpr PetscInt numAdditionalFields = 0;
  IncompressibleCascaded2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  IncompressibleCascaded2d(Real _om) : omega(_om) {}
  IncompressibleCascaded2d() : omega(1.) {}
  IncompressibleCascaded2d<D2Q9,CodiReverseType<Real>> getCodiAdjoint() const
  {
    return IncompressibleCascaded2d<D2Q9,CodiReverseType<Real>>(this->omega);
  }
  ReverseIncompressibleCascaded2d<D2Q9,Real> getSourceAdjoint() const
  {
    return ReverseIncompressibleCascaded2d<D2Q9,Real>(this->omega);
  }

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac) const
  {
    Real k[D2Q9::numDOF - 3];
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqr = mac[1]*mac[1] + mac[2]*mac[2];
    /* This is the fun part! */
    k[0] = (mac[0]*uSqr - fdistIn[8] - fdistIn[6] - fdistIn[2] - fdistIn[4] -
	    2.*(fdistIn[1] + fdistIn[3] + fdistIn[7] + fdistIn[5] - mac[0]/3.))/12.;
    k[1] = omega*(fdistIn[6] + fdistIn[2] - fdistIn[8] - fdistIn[4] + mac[0]*
		  (mac[1]*mac[1] - mac[2]*mac[2]))/4.;
    k[2] = omega*((fdistIn[7] + fdistIn[3] - fdistIn[5] - fdistIn[1]) -
		  mac[0]*mac[1]*mac[2])/4.;
    k[3] = (-((fdistIn[1] + fdistIn[3] - fdistIn[7] - fdistIn[5]
	       - 2.*mac[0]*mac[1]*mac[1]*mac[2] +
	       mac[2]*(mac[0] - fdistIn[6] - fdistIn[2] - fdistIn[0]))/4. +
	      0.5*mac[1]*(fdistIn[7] - fdistIn[5] - fdistIn[1] + fdistIn[3]))) +
      0.5*mac[2]*(-3.*k[0] - k[1]) + 2.*mac[1]*k[2];
    k[4] = (-((fdistIn[3] + fdistIn[5] - fdistIn[1] - fdistIn[7] -
	       2.*mac[0]*mac[2]*mac[2]*mac[1] +
	       mac[1]*(mac[0] - fdistIn[4] - fdistIn[8] - fdistIn[0]))/4. +
	      0.5*mac[2]*(fdistIn[7] + fdistIn[3] - fdistIn[1] - fdistIn[5]))) +
      0.5*mac[1]*(-3.*k[0] + k[1]) + 2.*mac[2]*k[2];
    k[5] = (0.25*(mac[0]/9. - fdistIn[7] - fdistIn[5] - fdistIn[1] - fdistIn[3] +
		  2.*(mac[1]*(fdistIn[7] - fdistIn[5] + fdistIn[1] - fdistIn[3]) +
		      mac[2]*(fdistIn[7] + fdistIn[5] - fdistIn[1] - fdistIn[3])) +
		  4.*mac[1]*mac[2]*(fdistIn[5] - fdistIn[7] + fdistIn[1] - fdistIn[3]) -
		  mac[1]*mac[1]*(fdistIn[6] + fdistIn[7] + fdistIn[5]
				 + fdistIn[2] + fdistIn[1] + fdistIn[3]) +
		  mac[2]*mac[2]*(3.*mac[0]*mac[1]*mac[1] - fdistIn[8] -
				 fdistIn[7] - fdistIn[5] - fdistIn[1] -
				 fdistIn[3] - fdistIn[4]))) -
      2.*k[0] - 2.*mac[1]*k[4] - 2.*mac[2]*k[3] + 4.*mac[1]*mac[2]*k[2] -
      1.5*k[0]*uSqr + 0.5*k[1]*(mac[1]*mac[1] - mac[2]*mac[2]);

    momentToDistribution(fdistIn,fdistOut,k);
  }
private:
  void momentToDistribution(const Real* fdistIn, Real* fdistOut,
                            PETSC_RESTRICT Real* k) const
  {
    fdistOut[0] = fdistIn[0] - 4.*k[0] + 4.*k[5];
    fdistOut[1] = fdistIn[1] + 2.*k[0] + k[2]
      + k[3] - k[4] + k[5];
    fdistOut[2] = fdistIn[2] - k[0] - k[1] - 2.*k[3] - 2.*k[5];
    fdistOut[3] = fdistIn[3] + 2.*k[0] - k[2]
      + k[3] + k[4] + k[5];
    fdistOut[4] = fdistIn[4] - k[0] + k[1] - 2.*k[4] - 2.*k[5];
    fdistOut[5] = fdistIn[5] + 2.*k[0] + k[2] - k[3]
      + k[4] + k[5];
    fdistOut[6] = fdistIn[6] - k[0] - k[1]
      + 2.*k[3] - 2.*k[5];
    fdistOut[7] = fdistIn[7]+ 2.*k[0]
      - k[2] - k[3] - k[4] + k[5];
    fdistOut[8] = fdistIn[8] - k[0] + k[1] + 2.*k[4] - 2.*k[5];
  }
  Real omega;
};

#endif
