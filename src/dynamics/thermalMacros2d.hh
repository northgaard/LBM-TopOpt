#ifndef THERMALMACROS2D
#define THERMALMACROS2D

#include "petsc.h"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
struct StandardThermalMacros2d {

  static void compute(const Real* fdist, Real& T){
    T = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      T += fdist[ii];
    }
  }
};

/* Lattice specializations */
template <class Real>
struct StandardThermalMacros2d<D2Q4,Real> {

  static void compute(const Real* fdist, Real& T){
    T = fdist[0] + fdist[1] + fdist[2] + fdist[3];
  }
};

template <class Real>
struct StandardThermalMacros2d<D2Q5,Real> {

  static void compute(const Real* fdist, Real& T){
    T = fdist[0] + fdist[1] + fdist[2] + fdist[3] + fdist[4];
  }
};

template <class Real>
struct StandardThermalMacros2d<D2Q5_TMRT,Real> {

  static void compute(const Real* fdist, Real& T){
    T = fdist[0] + fdist[1] + fdist[2] + fdist[3] + fdist[4];
  }
};

#endif
