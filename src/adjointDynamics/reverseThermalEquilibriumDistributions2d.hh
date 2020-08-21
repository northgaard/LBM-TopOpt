#ifndef REVERSETHERMALEQDIST2D
#define REVERSETHERMALEQDIST2D

#include "petsc.h"

template <class Lattice, class Real = PetscScalar>
struct ReverseThermalFeq2d {

  static void eq(PetscInt id, const Real ux, const Real uy, const Real T,
                 Real& uxRev, Real& uyRev, Real& TRev, const Real eqRev)
  {
    Real sc = Lattice::ex[id]*ux + Lattice::ey[id]*uy;
    Real temp = T*Lattice::weights[id]*Lattice::csSqInv*eqRev;
    uxRev += Lattice::ex[id]*temp;
    uyRev += Lattice::ey[id]*temp;
    TRev += Lattice::weights[id]*(1. + Lattice::csSqInv*sc)*eqRev;
  }
};

#endif
