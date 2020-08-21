#ifndef PBBCOLLISION
#define PBBCOLLISION

#include "petsc.h"
#include "core/codiheader.hh"
#include "core/meta.hh"
#include "adjointDynamics/reversePartialBounceback.hh"
#include <type_traits>

template <class BaseCollision, class InterpolationFunction>
class PartialBouncebackCollision {

  template <class T, class U>
  using EnableCodi = typename std::enable_if<(has_codi_adjoint<T>::value &&
                                              has_codi_adjoint<U>::value)>::type;
  template <class T, class U>
  using EnableSource = typename std::enable_if<(has_source_adjoint<T>::value &&
                                                has_diff_function<U>::value)>::type;
  BaseCollision baseOp;
  InterpolationFunction interp;
public:
  using LatticeType = typename BaseCollision::LatticeType;
  using RealType = typename BaseCollision::RealType;
  using Equilibrium = typename BaseCollision::Equilibrium;
  using Macros = typename BaseCollision::Macros;
  static constexpr PetscInt numAdditionalFields = 1;

  PartialBouncebackCollision(BaseCollision _op, InterpolationFunction _in)
    : baseOp(_op), interp(_in){}
  PartialBouncebackCollision() : baseOp(), interp(){}
  template <class T = BaseCollision, class U = InterpolationFunction>
  auto getCodiAdjoint(EnableCodi<T,U>* = nullptr) const
    -> PartialBouncebackCollision<decltype(baseOp.getCodiAdjoint()),
                                  decltype(interp.getCodiAdjoint())>
  {
    auto codiInter = interp.getCodiAdjoint();
    auto codiOp = baseOp.getCodiAdjoint();
    return PartialBouncebackCollision<decltype(codiOp),decltype(codiInter)>
      (codiOp,codiInter);
  }
  template <class T = BaseCollision, class U = InterpolationFunction>
  auto getSourceAdjoint(EnableSource<T,U>* = nullptr) const
    -> ReversePartialBouncebackCollision<BaseCollision,
                                         decltype(baseOp.getSourceAdjoint()),
                                         InterpolationFunction>
  {
    auto reverseBaseOp = baseOp.getSourceAdjoint();
    return ReversePartialBouncebackCollision<BaseCollision,decltype(reverseBaseOp),
                                             InterpolationFunction>
      (baseOp,reverseBaseOp,interp);
  }

  inline void operator()(const RealType* fdistIn, RealType* fdistOut,
                         RealType* mac, RealType* obst)
  {
    RealType ftemp, imult;
    baseOp(fdistIn,fdistOut,mac);
    imult = 0.5*interp(*obst);

    for (size_t dd = 1; dd <= LatticeType::half; ++dd){
      ftemp = fdistOut[dd];
      fdistOut[dd] += imult*(fdistOut[dd+LatticeType::half] - fdistOut[dd]);
      fdistOut[dd+LatticeType::half] += imult*(ftemp - fdistOut[dd+LatticeType::half]);
    }
  }
};

#endif
