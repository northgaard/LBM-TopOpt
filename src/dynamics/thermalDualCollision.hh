#ifndef THERMALDUALCOLLISION
#define THERMALDUALCOLLISION

#include "petsc.h"
#include "core/codiheader.hh"
#include "adjointDynamics/reverseThermalDualCollision.hh"
#include <type_traits>

template <class IsoCollision, class ThermalCollision>
class ThermalDualCollision {

  // Assume shared additional fields for now
  static_assert(IsoCollision::numAdditionalFields ==
                ThermalCollision::numAdditionalFields,
                "Number of additional fields must be the same.\n");
  static_assert(std::is_same<typename IsoCollision::RealType,
                typename ThermalCollision::RealType>::value,
                "Floating point type of collision operators must be the same.\n");
  using IsoCollisionType = IsoCollision;
  using ThermalCollisionType = ThermalCollision;
  IsoCollision isoOp;
  ThermalCollision thermalOp;
public:
  using IsoLatticeType = typename IsoCollision::LatticeType;
  using ThermalLatticeType = typename ThermalCollision::LatticeType;
  using RealType = typename IsoCollision::RealType;
  using IsoEquilibriumType = typename IsoCollision::Equilibrium;
  using ThermalEquilibriumType = typename ThermalCollision::Equilibrium;
  using IsoMacrosType = typename IsoCollision::Macros;
  using ThermalMacrosType = typename ThermalCollision::Macros;
  static constexpr PetscInt numAdditionalFields = IsoCollision::numAdditionalFields;

  ThermalDualCollision(IsoCollision _iop, ThermalCollision _top)
    : isoOp(_iop), thermalOp(_top){}
  ThermalDualCollision() : isoOp(), thermalOp(){}
  auto getCodiAdjoint() const
    -> ThermalDualCollision<decltype(isoOp.getCodiAdjoint()),
                            decltype(thermalOp.getCodiAdjoint())>
  {
    auto codiIsoOp = isoOp.getCodiAdjoint();
    auto codiThermalOp = thermalOp.getCodiAdjoint();
    return ThermalDualCollision<decltype(codiIsoOp),decltype(codiThermalOp)>
      (codiIsoOp,codiThermalOp);
  }
  auto getSourceAdjoint() const
    -> ReverseThermalDualCollision<IsoCollision,ThermalCollision,
                                   decltype(isoOp.getSourceAdjoint()),
                                   decltype(thermalOp.getSourceAdjoint())>
  {
    auto reverseIsoOp = isoOp.getSourceAdjoint();
    auto reverseThermalOp = thermalOp.getSourceAdjoint();
    return ReverseThermalDualCollision<IsoCollision,ThermalCollision,
                                       decltype(reverseIsoOp),decltype(reverseThermalOp)>
      (reverseIsoOp,reverseThermalOp);
  }

  void operator()(const RealType* fdistIn, RealType* fdistOut,
                         RealType* mac, RealType* obst)
  {
    isoOp(fdistIn,fdistOut,mac,obst);
    thermalOp(fdistIn + IsoLatticeType::numDOF,
              fdistOut + IsoLatticeType::numDOF,
              mac, obst);
  }
};

#endif
