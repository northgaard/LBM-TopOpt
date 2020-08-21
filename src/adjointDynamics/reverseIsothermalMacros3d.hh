#ifndef REVERSEISOTHERMALMACROS3D
#define REVERSEISOTHERMALMACROS3D

#include "petsc.h"

template <class Lattice, class Real = PetscScalar>
struct ReverseIncompressibleMacros3d {

  static void compute(Real* fdistRev, const Real rhoRev, const Real uxRev,
                      const Real uyRev, const Real uzRev)
  {
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistRev[ii] += rhoRev + Lattice::ex[ii]*uxRev +
        Lattice::ey[ii]*uyRev + Lattice::ez[ii]*uzRev;
    }
  }
};

#endif
