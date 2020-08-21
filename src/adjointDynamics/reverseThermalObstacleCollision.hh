#ifndef REVERSETHERMALOBSTACLECOLLISION
#define REVERSETHERMALOBSTACLECOLLISION

#include "petsc.h"
#include <type_traits>

template <class BaseOperator, class ReverseBaseOperator,
          class InterpolationFunction>
class ReverseThermalObstacleCollision {

  static_assert(std::is_same<typename BaseOperator::RealType,
                typename ReverseBaseOperator::RealType>::value,
                "Base and reverse operator must have consistent real type.\n");
public:
  using LatticeType = typename BaseOperator::LatticeType;
  using RealType = typename BaseOperator::RealType;

  ReverseThermalObstacleCollision(const BaseOperator& _bop, ReverseBaseOperator _rop,
                                  InterpolationFunction _in)
    : reverseBaseOp(_rop), interp(_in)
  {
    refDiffusivity = _bop.getDiffusivity();
  }

  inline void operator()(RealType* fdistIn, RealType* fdistOut, RealType* mac,
                         RealType* obst, RealType* fdistInRev, RealType* fdistOutRev,
                         RealType* macRev, RealType* obstRev)
  {
    RealType imult = refDiffusivity * interp(*obst);
    RealType temp = LatticeType::csSqInv*imult + 0.5;
    RealType omega = 1./temp;
    RealType omegaRev = 0.;
    reverseBaseOp(fdistIn,fdistOut,mac,fdistInRev,fdistOutRev,macRev,omega,&omegaRev);
    RealType imultRev = -(LatticeType::csSqInv*omegaRev/(temp*temp));
    *obstRev += refDiffusivity*interp.diff(*obst)*imultRev;
  }
private:
  ReverseBaseOperator reverseBaseOp;
  InterpolationFunction interp;
  RealType refDiffusivity;
};

#endif
