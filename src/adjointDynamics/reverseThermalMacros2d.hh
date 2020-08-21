#ifndef REVERSETHERMALMACROS2D
#define REVERSETHERMALMACROS2D

#include "petsc.h"

template <class Lattice, class Real = PetscScalar>
struct ReverseThermalMacros2d {

  static void compute(Real* fdistRev, const Real TRev)
  {
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdistRev[ii] += TRev;
    }
  }
};

#endif
