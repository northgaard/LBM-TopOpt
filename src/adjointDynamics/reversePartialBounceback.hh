#ifndef REVERSEPBBCOLLISION
#define REVERSEPBBCOLLISION

#include "petsc.h"

template <class BaseOperator, class ReverseBaseOperator,
          class InterpolationFunction>
class ReversePartialBouncebackCollision {

public:
  using Lattice = typename BaseOperator::LatticeType;
  using RealType = typename BaseOperator::RealType;

  ReversePartialBouncebackCollision(BaseOperator _op, ReverseBaseOperator _rop,
                                    InterpolationFunction _in)
    : baseOp(_op), reverseBaseOp(_rop), interp(_in){}

  inline void operator()(RealType* fdistIn, RealType* fdistOut, RealType* mac,
                         RealType* obst, RealType* fdistInRev, RealType* fdistOutRev,
                         RealType* macRev, RealType* obstRev)
  {
    RealType imult;
    RealType ftempRev;
    RealType imultRev = 0.;

    baseOp(fdistIn,fdistOut,mac);
    imult = 0.5*interp(*obst);
    for (size_t dd = 1; dd <= Lattice::half; ++dd){
      imultRev += (fdistOut[dd] - fdistOut[dd+Lattice::half]) *
        fdistOutRev[dd+Lattice::half];
      imultRev += (fdistOut[dd+Lattice::half] - fdistOut[dd]) *
        fdistOutRev[dd];
      ftempRev = fdistOutRev[dd+Lattice::half];
      fdistOutRev[dd+Lattice::half] +=
        imult*(fdistOutRev[dd] - fdistOutRev[dd+Lattice::half]);
      fdistOutRev[dd] += imult*(ftempRev - fdistOutRev[dd]);
    }

    *obstRev += 0.5*interp.diff(*obst)*imultRev;
    reverseBaseOp(fdistIn,fdistOut,mac,fdistInRev,fdistOutRev,macRev);
  }
private:
  BaseOperator baseOp;
  ReverseBaseOperator reverseBaseOp;
  InterpolationFunction interp;
};

#endif
