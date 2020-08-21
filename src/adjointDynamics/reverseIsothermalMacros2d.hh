#ifndef REVERSEISOTHERMALMACROS2D
#define REVERSEISOTHERMALMACROS2D

#include "petsc.h"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
struct ReverseIncompressibleMacros2d {

  static void compute(Real* fdistRev, const Real rhoRev, const Real uxRev,
		      const Real uyRev)
  {
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistRev[ii] += rhoRev + Lattice::ex[ii]*uxRev +
        Lattice::ey[ii]*uyRev;
    }
  }
};

/* D2Q9 optimization */
template <class Real>
struct ReverseIncompressibleMacros2d<D2Q9,Real> {

  static void compute(Real* fdistRev, const Real rhoRev, const Real uxRev,
                      const Real uyRev)
  {
    fdistRev[0] += rhoRev;
    fdistRev[1] += rhoRev + uxRev - uyRev;
    fdistRev[2] += rhoRev - uyRev;
    fdistRev[3] += rhoRev - uxRev - uyRev;
    fdistRev[4] += rhoRev - uxRev;
    fdistRev[5] += rhoRev - uxRev + uyRev;
    fdistRev[6] += rhoRev + uyRev;
    fdistRev[7] += rhoRev + uxRev + uyRev;
    fdistRev[8] += rhoRev + uxRev;
  }
};

#endif
