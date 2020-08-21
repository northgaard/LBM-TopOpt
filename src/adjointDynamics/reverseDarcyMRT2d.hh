#ifndef REVERSEDARCYMRT2D
#define REVERSEDARCYMRT2D

#include "adjointDynamics/reverseIncompressibleMRT2d.hh"

template <class Lattice, class InterpolationFunction, class Real = PetscScalar>
class ReverseDarcyMRT2d {};

template <class InterpolationFunction, class Real>
class ReverseDarcyMRT2d<D2Q9,InterpolationFunction,Real> {

public:
  using Lattice = D2Q9;
  using RevMacros = ReverseIncompressibleMacros2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  using RealType = Real;

  ReverseDarcyMRT2d(ReverseIncompressibleMRT2d<D2Q9,Real> _revOp,
                    InterpolationFunction _in, Real _om, Real _fx,
                    Real _fy, Real _am)
    : reverseBaseOp(_revOp), interp(_in), omega(_om), cfx(_fx), cfy(_fy), amax(_am)
  {}

  void operator()(const Real* fdistIn, Real* fdistOut, Real* mac,
                  Real* obst, Real* fdistInRev, Real* fdistOutRev,
                  Real* macRev, Real* obstRev)
  {

  }

private:
  static constexpr PetscScalar s1 = 1.4;
  static constexpr PetscScalar s2 = 1.4;
  static constexpr PetscScalar s46 = 1.2;
  ReverseIncompressibleMRT2d<D2Q9,Real> reverseBaseOp;
  InterpolationFunction interp;
  Real omega;
  const Real cfx, cfy, amax;
  // Fractions
  static constexpr PetscScalar oo3 = 1./3.;
  static constexpr PetscScalar oo6 = 1./6.;
  static constexpr PetscScalar oo9 = 1./9.;
  static constexpr PetscScalar oo12 = 1./12.;
  static constexpr PetscScalar oo18 = 1./18.;
  static constexpr PetscScalar oo24 = 1./24.;
  static constexpr PetscScalar oo36 = 1./36.;
};

#endif
