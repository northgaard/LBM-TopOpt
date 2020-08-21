#ifndef THERMALPBB2D
#define THERMALPBB2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/codiheader.hh"
#include <type_traits>

template <class IsoLattice, class ThermalLattice,
	  template <class> class IsoInterpolation,
	  template <class> class ThermalInterpolation,
	  template <class,class> class IsoCollision,
	  template <class,class> class ThermalCollision,
	  class Real = PetscScalar>
class ThermalPartialBounceback2d {

  template <class T>
  using Partial = ThermalPartialBounceback2d
    <IsoLattice,ThermalLattice,IsoInterpolation,
     ThermalInterpolation,IsoCollision,ThermalCollision,T>;

  template <class,class,
	    template <class> class,
	    template <class> class,
	    template <class,class> class,
	    template <class,class> class,
	    class OtherReal>
  friend class ThermalPartialBounceback2d;

public:

  ThermalPartialBounceback2d(ThermalFlowParameters par,
			   IsoInterpolation<Real> _inter,
			   ThermalInterpolation<Real> _tinter)
    : colOp(par.incPar), thermOp(par),
      interp(_inter), thermInterp(_tinter)
  {
    refDiffusivity = (par.incPar.velocityChar * par.incPar.lengthChar) /
      (par.incPar.ReynoldsNumber * par.PrandtlNumber);
  }

  ThermalPartialBounceback2d()
    : colOp(), thermOp(), interp(), thermInterp() {}

  template <class OtherReal>
  ThermalPartialBounceback2d(const Partial<OtherReal>& rhs)
    : colOp(rhs.colOp), thermOp(rhs.thermOp),
      interp(rhs.interp), thermInterp(rhs.thermInterp)
  {
    refDiffusivity = rhs.refDiffusivity;
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
		       Real* mac, Real obst) const
  {

    Real ftemp,imult;
    colOp(fdistIn,fdistOut,mac);

    imult = 0.5*interp(obst);

    for (size_t dd = 1; dd <= IsoLattice::half; ++dd){
      ftemp = fdistOut[dd];
      fdistOut[dd] +=
	imult*(fdistOut[dd+IsoLattice::half] - fdistOut[dd]);
      fdistOut[dd+IsoLattice::half] +=
	imult*(ftemp - fdistOut[dd+IsoLattice::half]);
    }

    imult = refDiffusivity * thermInterp(obst);
    thermOp.setOmega(1./(2.*imult + 0.5));
    thermOp(fdistIn + IsoLattice::numDOF,
	    fdistOut + IsoLattice::numDOF,
	    mac);
    
  }

private:

  IsoCollision<IsoLattice,Real> colOp;
  ThermalCollision<ThermalLattice,Real> thermOp;
  IsoInterpolation<Real> interp;
  ThermalInterpolation<Real> thermInterp;
  PetscScalar refDiffusivity;

};

template <class IsoLattice, class ThermalLattice,
	  template <class> class IsoInterpolation,
	  template <class> class ThermalInterpolation,
	  template <class,class> class IsoCollision,
	  template <class,class> class ThermalCollision,
	  class Real = PetscScalar>
class HeatingThermalPartialBounceback2d {

  template <class T>
  using Partial = HeatingThermalPartialBounceback2d
    <IsoLattice,ThermalLattice,IsoInterpolation,
     ThermalInterpolation,IsoCollision,ThermalCollision,T>;

  template <class,class,
	    template <class> class,
	    template <class> class,
	    template <class,class> class,
	    template <class,class> class,
	    class OtherReal>
  friend class HeatingThermalPartialBounceback2d;

public:

  HeatingThermalPartialBounceback2d(ThermalFlowParameters par,
				    PetscScalar _hv,
				    IsoInterpolation<Real> _inter,
				    ThermalInterpolation<Real> _tinter)
    : colOp(par.incPar), thermOp(par),
      interp(_inter), thermInterp(_tinter), heatingValue(_hv)
  {
    refDiffusivity = (par.incPar.velocityChar * par.incPar.lengthChar) /
      (par.incPar.ReynoldsNumber * par.PrandtlNumber);
  }

  HeatingThermalPartialBounceback2d()
    : colOp(), thermOp(), interp(), thermInterp() {}

  template <class OtherReal>
  HeatingThermalPartialBounceback2d(const Partial<OtherReal>& rhs)
    : colOp(rhs.colOp), thermOp(rhs.thermOp),
      interp(rhs.interp), thermInterp(rhs.thermInterp)
  {
    heatingValue = rhs.heatingValue;
    refDiffusivity = rhs.refDiffusivity;
  }

  inline void operator()(const Real* fdistIn, Real* fdistOut,
		       Real* mac, Real obst) const
  {

    Real ftemp,imult;
    colOp(fdistIn,fdistOut,mac);

    imult = 0.5*interp(obst);

    for (size_t dd = 1; dd <= IsoLattice::half; ++dd){
      ftemp = fdistOut[dd];
      fdistOut[dd] +=
	imult*(fdistOut[dd+IsoLattice::half] - fdistOut[dd]);
      fdistOut[dd+IsoLattice::half] +=
	imult*(ftemp - fdistOut[dd+IsoLattice::half]);
    }

    imult = refDiffusivity * thermInterp(obst);
    thermOp.setOmega(1./(2.*imult + 0.5));
    thermOp(fdistIn + IsoLattice::numDOF,
	    fdistOut + IsoLattice::numDOF,
	    mac);

    // Obstacle heating
    for (size_t dd = 0; dd < ThermalLattice::numDOF; ++dd){
      fdistOut[dd + IsoLattice::numDOF] +=
	(1. - obst)*ThermalLattice::weights[dd]*heatingValue;
    }
  }

private:

  IsoCollision<IsoLattice,Real> colOp;
  ThermalCollision<ThermalLattice,Real> thermOp;
  IsoInterpolation<Real> interp;
  ThermalInterpolation<Real> thermInterp;
  PetscScalar heatingValue;
  PetscScalar refDiffusivity;

public:

    using CodiReverseType = typename std::conditional<
    std::is_same<Real,PetscScalar>::value,
    Partial<ReversePetscScalar>,
    void>::type;

};

#endif
