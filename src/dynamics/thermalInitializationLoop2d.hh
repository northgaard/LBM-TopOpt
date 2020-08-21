#ifndef THERMALINITLOOP2D
#define THERMALINITLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"

template <class IsothermalLattice,
	  class IsothermalEq,
	  class ThermalEq>
class ThermalEquilibriumInitializationLoop2d {

public:

  ThermalEquilibriumInitializationLoop2d(Box2d _box) : boundingBox(_box) {}

  void operator()(PetscScalar*** fdist, PetscScalar*** mac) const
  {

    PetscScalar uSqr;
    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
	uSqr = mac[jj][ii][1]*mac[jj][ii][1] + mac[jj][ii][2]*mac[jj][ii][2];
	IsothermalEq::setAllEquilibria(fdist[jj][ii],mac[jj][ii][0],
				       mac[jj][ii][1],mac[jj][ii][2],uSqr);
	ThermalEq::setAllEquilibria(fdist[jj][ii] + IsothermalLattice::numDOF,
				  mac[jj][ii][1],mac[jj][ii][2],mac[jj][ii][3]);
      }
    }
  }

private:

  Box2d boundingBox;

};

#endif
