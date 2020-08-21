#ifndef THERMALEQDIST2D
#define THERMALEQDIST2D

#include "petsc.h"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
struct ThermalFeq2d {

  static Real eq(PetscInt id, const Real ux, const Real uy,
		 const Real T)
  {
    Real sc = Lattice::ex[id]*ux + Lattice::ey[id]*uy;
    return Lattice::weights[id]*T*(1. + Lattice::csSqInv*sc);
  }

  static void setAllEquilibria(Real* fdist, const Real ux,
			       const Real uy, const Real T)
  {
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdist[ii] = ThermalFeq2d<Lattice,Real>::eq(ii,ux,uy,T);
    }
  }
};

template <class Lattice, class Real = PetscScalar>
struct ThermalMRTFeq2d {};

template <class Real>
struct ThermalMRTFeq2d<D2Q5_TMRT,Real> {

  static void setAllEquilibria(Real* fdist, const Real ux,
                               const Real uy, const Real T)
  {
    Real moments[D2Q5_TMRT::numDOF - 1];
    moments[0] = ux*T;
    moments[1] = uy*T;
    moments[2] = -2.*T;
    moments[3] = 0.;
    // Map equilibrium moments to distribution space
    fdist[0] = 0.2*T - 0.2*moments[2];
    fdist[1] = 0.2*T - 0.5*moments[1] + 0.05*moments[2] - 0.25*moments[3];
    fdist[2] = 0.2*T - 0.5*moments[0] + 0.05*moments[2] + 0.25*moments[3];
    fdist[3] = 0.2*T + 0.5*moments[1] + 0.05*moments[2] - 0.25*moments[3];
    fdist[4] = 0.2*T + 0.5*moments[0] + 0.05*moments[2] + 0.25*moments[3];
  }
};

#endif
