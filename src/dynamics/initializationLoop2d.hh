#ifndef INITIALIZATIONLOOP2D
#define INITIALIZATIONLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"

template <class Equilibrium>
class EquilibriumInitializationLoop2d {

public:

  EquilibriumInitializationLoop2d(Box2d);
  void operator()(void*,void*,void*) const;

private:

  Box2d boundingBox;

};

template <class Equilibrium>
EquilibriumInitializationLoop2d<Equilibrium>
::EquilibriumInitializationLoop2d(Box2d _box) :  boundingBox(_box) {}

template <class Equilibrium>
void EquilibriumInitializationLoop2d<Equilibrium>::
operator()(void* fdist_v, void* mac_v, void*) const
{
  
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscScalar*** mac = (PetscScalar***) mac_v;
  PetscScalar uSqr;

  for (auto jj : boundingBox.yRange){
    for (auto ii : boundingBox.xRange){
      uSqr = mac[jj][ii][1]*mac[jj][ii][1] + mac[jj][ii][2]*mac[jj][ii][2];
      Equilibrium::setAllEquilibria(fdist[jj][ii],mac[jj][ii][0],mac[jj][ii][1],
				    mac[jj][ii][2],uSqr);
    }
  }
}

#endif
