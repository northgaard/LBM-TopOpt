#ifndef THERMALOBSTACLECOLLISION
#define THERMALOBSTACLECOLLISION

#include "petsc.h"
#include "core/codiheader.hh"
#include "adjointDynamics/reverseThermalObstacleCollision.hh"

template <class BaseCollision, class InterpolationFunction>
class ThermalObstacleCollision {

public:
  using LatticeType = typename BaseCollision::LatticeType;
  using RealType = typename BaseCollision::RealType;
  using Equilibrium = typename BaseCollision::Equilibrium;
  using Macros = typename BaseCollision::Macros;
  static constexpr PetscInt numAdditionalFields = 1;
private:
  BaseCollision baseOp;
  InterpolationFunction interp;
  RealType refDiffusivity;
public:
  ThermalObstacleCollision(BaseCollision _op, InterpolationFunction _in)
    : baseOp(_op), interp(_in)
  {
    refDiffusivity = baseOp.getDiffusivity();
  }
  ThermalObstacleCollision() : baseOp(), interp(){}
  auto getCodiAdjoint() const
    -> ThermalObstacleCollision<decltype(baseOp.getCodiAdjoint()),
                                decltype(interp.getCodiAdjoint())>
  {
    auto codiInter = interp.getCodiAdjoint();
    auto codiOp = baseOp.getCodiAdjoint();
    return ThermalObstacleCollision<decltype(codiOp),decltype(codiInter)>
      (codiOp,codiInter);
  }
  auto getSourceAdjoint() const
    -> ReverseThermalObstacleCollision<BaseCollision,decltype(baseOp.getSourceAdjoint()),
                                       InterpolationFunction>
  {
    auto reverseBaseOp = baseOp.getSourceAdjoint();
    return ReverseThermalObstacleCollision<BaseCollision,decltype(reverseBaseOp),
                                           InterpolationFunction>
      (baseOp,reverseBaseOp,interp);
  }

  inline void operator()(const RealType* fdistIn, RealType* fdistOut,
                         RealType* mac, RealType* obst)
  {
    RealType imult = refDiffusivity * interp(*obst);
    RealType omega = 1./(LatticeType::csSqInv*imult + 0.5);
    baseOp(fdistIn,fdistOut,mac,omega);
  }
};

#endif
