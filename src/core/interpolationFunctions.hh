#ifndef INTERPOLATIONFUNC
#define INTERPOLATIONFUNC

#include "petsc.h"
#include "core/codiheader.hh"

template <class Real = PetscScalar>
class BPInterpolation {

public:
  BPInterpolation() : penal(1.) {}
  BPInterpolation(PetscScalar _p) : penal(_p) {}
  BPInterpolation<CodiReverseType<Real>> getCodiAdjoint() const
  {
    return BPInterpolation<CodiReverseType<Real>>(penal);
  }

  Real operator()(const Real s) const
  {
    return 1. - s*(1. + penal)/(penal + s);
  }
  Real diff(const Real s) const
  {
    return (-(1. + penal)/(penal + s) + s*(1. + penal)/((penal + s)*(penal + s)));
  }
private:
  PetscScalar penal;
};

template <class Real = PetscScalar>
class RampInterpolation {

public:
  RampInterpolation() : penal(1.), Ck(1.) {}
  RampInterpolation(PetscScalar _p, PetscScalar _C)
    : penal(_p), Ck(_C) {}
  RampInterpolation<CodiReverseType<Real>> getCodiAdjoint() const
  {
    return RampInterpolation<CodiReverseType<Real>>(penal,Ck);
  }

  Real operator()(const Real s) const
  {
    return (s*(Ck*(1. + penal) - 1.) + 1.) / (Ck*(1. + penal*s));
  }
  Real diff(const Real s) const
  {
    return (Ck*(1. + penal) - 1.) / (Ck*(1. + penal*s))
      - ((s*(Ck*(1. + penal) - 1.) + 1.)*penal) /
      (Ck*(1. + penal*s)*(1. + penal*s));
  }
private:
  PetscScalar penal,Ck;
};

#endif
