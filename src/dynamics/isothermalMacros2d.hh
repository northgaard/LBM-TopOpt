#ifndef ISOTHERMALMACROS2D
#define ISOTHERMALMACROS2D

#include "petsc.h"
#include "latticeBoltzmann/lattices2d.hh"

template <class Lattice, class Real = PetscScalar>
struct IncompressibleMacros2d {

  static void compute(const Real* fdist, Real& rho, Real& ux, Real& uy)
  {
    rho = 0.;
    ux = 0.;
    uy = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
      ux += Lattice::ex[ii]*fdist[ii];
      uy += Lattice::ey[ii]*fdist[ii];
    }
  }
  static void computeRho(const Real* fdist, Real& rho)
  {
    rho = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
    }
  }
  static void computeUx(const Real* fdist, Real& ux)
  {
    ux = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      ux += Lattice::ex[ii]*fdist[ii];
    }
  }
  static void computeUy(const Real* fdist, Real& uy)
  {
    uy = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      uy += Lattice::ey[ii]*fdist[ii];
    }
  }
};

template <class Real>
struct IncompressibleMacros2d<D2Q9,Real> {

  static void compute(const Real* fdist, Real& rho, Real& ux, Real& uy){
    ux = fdist[1] + fdist[7] + fdist[8];
    uy = fdist[5] + fdist[6] + fdist[7];
    rho = ux + uy - fdist[7] + fdist[0] + fdist[2] + fdist[3] + fdist[4];
    ux -= fdist[3] + fdist[4] + fdist[5];
    uy -= fdist[1] + fdist[2] + fdist[3];
  }
  static void computeRho(const Real* fdist, Real& rho)
  {
    rho = fdist[0]+fdist[1]+fdist[2]+fdist[3]+fdist[4]+fdist[5]+fdist[6]
      +fdist[7]+fdist[8];
  }
  static void computeUx(const Real* fdist, Real& ux)
  {
    ux = fdist[1]+fdist[7]+fdist[8]-(fdist[3]+fdist[4]+fdist[5]);
  }
  static void computeUy(const Real* fdist, Real& uy)
  {
    uy = fdist[5]+fdist[6]+fdist[7]-(fdist[1]+fdist[2]+fdist[3]);
  }
};

template <class Lattice, class Real = PetscScalar>
struct StandardMacros2d {

  static void compute(const Real* fdist, Real& rho, Real& ux, Real& uy){
    rho = 0.;
    ux = 0.;
    uy = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
      ux += Lattice::ex[ii]*fdist[ii];
      uy += Lattice::ey[ii]*fdist[ii];
    }
    Real inv = 1./rho;
    ux *= inv; uy *= inv;
  }
  static void computeRho(const Real* fdist, Real& rho)
  {
    rho = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      rho += fdist[ii];
    }
  }
  static void computeUx(const Real* fdist, Real& ux)
  {
    ux = 0.;
    PetscScalar rho = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      ux += Lattice::ex[ii]*fdist[ii];
      rho += fdist[ii];
    }
    ux /= rho;
  }
  static void computeUy(const Real* fdist, Real& uy)
  {
    uy = 0.;
    PetscScalar rho = 0.;
    for (size_t ii = 0; ii < Lattice::numDOF; ++ii){
      uy += Lattice::ey[ii]*fdist[ii];
      rho += fdist[ii];
    }
    uy /= rho;
  }
};

template <class Real>
struct StandardMacros2d<D2Q9,Real> {

  static void compute(const Real* fdist, Real& rho, Real& ux, Real& uy){
    ux = fdist[1] + fdist[7] + fdist[8];
    uy = fdist[5] + fdist[6] + fdist[7];
    rho = ux + uy - fdist[7] + fdist[0] + fdist[2] + fdist[3] + fdist[4];
    ux -= fdist[3] + fdist[4] + fdist[5];
    uy -= fdist[1] + fdist[2] + fdist[3];
    Real inv = 1./rho;
    ux *= inv; uy *= inv;
  }
  static void computeRho(const Real* fdist, Real& rho)
  {
    rho = fdist[0]+fdist[1]+fdist[2]+fdist[3]+fdist[4]+fdist[5]+fdist[6]
      +fdist[7]+fdist[8];
  }
  static void computeUx(const Real* fdist, Real& ux)
  {
    ux = fdist[1]+fdist[7]+fdist[8];
    PetscScalar rho = ux+fdist[0]+fdist[2]+fdist[3]+fdist[4]
      +fdist[5]+fdist[6];
    ux -= fdist[3]+fdist[4]+fdist[5];
    ux /= rho;
  }
  static void computeUy(const Real* fdist, Real& uy)
  {
    uy = fdist[5]+fdist[6]+fdist[7];
    PetscScalar rho = uy+fdist[0]+fdist[1]+fdist[2]+fdist[3]
      +fdist[4]+fdist[8];
    uy -= fdist[1]+fdist[2]+fdist[3];
    uy /= rho;
  }
};

#endif
