#ifndef ISOTHERMALMACROS3D
#define ISOTHERMALMACROS3D

#include "petsc.h"

template <class Lattice, class Real = PetscScalar>
struct IncompressibleMacros3d {

  static void compute(const Real* fdist, Real& rho,
                      Real& ux, Real& uy, Real& uz)
  {
    rho = 0.;
    ux = 0.;
    uy = 0.;
    uz = 0.;

    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
      ux += Lattice::ex[ii]*fdist[ii];
      uy += Lattice::ey[ii]*fdist[ii];
      uz += Lattice::ez[ii]*fdist[ii];
    }
  }
};

template <class Lattice, class Real = PetscScalar>
struct StandardMacros3d {

  static void compute(const Real* fdist, Real& rho,
                      Real& ux, Real& uy, Real& uz)
  {
    rho = 0.;
    ux = 0.;
    uy = 0.;
    uz = 0.;

    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
      ux += Lattice::ex[ii]*fdist[ii];
      uy += Lattice::ey[ii]*fdist[ii];
      uz += Lattice::ez[ii]*fdist[ii];
    }
    Real inv = 1./rho;
    ux *= inv; uy *= inv; uz *= inv;
  }
  static void computeMoments(const Real* fdist, Real& rho,
                             Real& jx, Real& jy, Real& jz)
  {
    rho = 0.;
    jx = 0.;
    jy = 0.;
    jz = 0.;

    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
      jx += Lattice::ex[ii]*fdist[ii];
      jy += Lattice::ey[ii]*fdist[ii];
      jz += Lattice::ez[ii]*fdist[ii];
    }
  }
};

#endif
