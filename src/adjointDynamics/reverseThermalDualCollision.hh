#ifndef REVERSETHERMALDUALCOLLISION
#define REVERSETHERMALDUALCOLLISION

#include "petsc.h"
#include "core/codiheader.hh"
#include <type_traits>

template <class IsoCollision, class ThermalCollision,
          class ReverseIsoCollision, class ReverseThermalCollision>
class ReverseThermalDualCollision {

public:
  using IsoLatticeType = typename IsoCollision::LatticeType;
  using ThermalLatticeType = typename ThermalCollision::LatticeType;
  using IsoMacrosType = typename IsoCollision::Macros;
  using ThermalMacrosType = typename ThermalCollision::Macros;
  using RealType = typename IsoCollision::RealType;

  ReverseThermalDualCollision(ReverseIsoCollision _ri, ReverseThermalCollision _rt)
    : rIsoOp(_ri), rThermalOp(_rt){}

  void operator()(RealType* fdistIn, RealType* fdistOut, RealType* mac,
                  RealType* obst, RealType* fdistInRev, RealType* fdistOutRev,
                  RealType* macRev, RealType* obstRev)
  {
    static constexpr PetscInt offset = IsoLatticeType::numDOF;
    IsoMacrosType::compute(fdistIn,mac[0],mac[1],mac[2]);
    rThermalOp(fdistIn + offset, fdistOut + offset,
               mac, obst, fdistInRev + offset, fdistOutRev + offset,
               macRev, obstRev);
    rIsoOp(fdistIn,fdistOut,mac,obst,
           fdistInRev,fdistOutRev,macRev,obstRev);
  }
private:
  ReverseIsoCollision rIsoOp;
  ReverseThermalCollision rThermalOp;
};

#endif
