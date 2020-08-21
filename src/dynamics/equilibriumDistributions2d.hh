#ifndef EQDIST2D
#define EQDIST2D

#include "petsc.h"

template <class Lattice, class Real = PetscScalar>
struct StandardFeq2d {

  static Real eq(PetscInt id, const Real rho, const Real ux,
		 const Real uy, const Real uSqr)
  {
    Real sc = Lattice::ex[id]*ux + Lattice::ey[id]*uy;
    return Lattice::weights[id]*rho*(1. + 3.*sc + 4.5*sc*sc - 1.5*uSqr);
  }
  static void setAllEquilibria(Real* fdist, const Real rho, const Real ux,
			       const Real uy, const Real uSqr)
  {
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdist[ii] = StandardFeq2d<Lattice,Real>::eq(ii,rho,ux,uy,uSqr);
    }
  }
};

template <class Lattice, class Real = PetscScalar>
struct IncompressibleFeq2d {

  static Real eq(PetscInt id, const Real rho, const Real ux,
			 const Real uy, const Real uSqr)
  {
    Real sc = Lattice::ex[id]*ux + Lattice::ey[id]*uy;
    return Lattice::weights[id]*(rho + 3.*sc + 4.5*sc*sc - 1.5*uSqr);
  }
  static void setAllEquilibria(Real* fdist, const Real rho, const Real ux,
			       const Real uy, const Real uSqr)
  {
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      fdist[ii] = IncompressibleFeq2d<Lattice,Real>::eq(ii,rho,ux,uy,uSqr);
    }
  }
};

#endif
