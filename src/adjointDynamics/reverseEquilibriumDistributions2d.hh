#ifndef REVERSEEQDIST2D
#define REVERSEEQDIST2D

#include "petsc.h"

template <class Lattice, class Real = PetscScalar>
struct ReverseIncompressibleFeq2d {

  static void eq(PetscInt id, const Real rho, const Real ux, const Real uy,
		 const Real uSqr, Real& rhoRev, Real& uxRev, Real& uyRev,
		 Real& uSqrRev, const Real eqRev)
  {
    Real sc = Lattice::ex[id]*ux + Lattice::ey[id]*uy;
    Real temp = Lattice::weights[id]*eqRev;
    Real scRev = (9.*sc + 3.)*temp;
    rhoRev += temp;
    uxRev += Lattice::ex[id]*scRev;
    uyRev += Lattice::ey[id]*scRev;
    uSqrRev -= 1.5*temp;
  }
};

#endif
