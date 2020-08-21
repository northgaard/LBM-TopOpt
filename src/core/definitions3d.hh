#ifndef DEFINITIONS3D
#define DEFINITIONS3D

#include "petsc.h"

template <class Real>
struct t_IsothermalMacros3d {
  Real rho;
  Real ux;
  Real uy;
  Real uz;
};

using IsothermalMacros3d = t_IsothermalMacros3d<PetscScalar>;

template <class Real>
struct t_ThermalMacros3d {
  Real rho;
  Real ux;
  Real uy;
  Real uz;
  Real T;
};

using ThermalMacros3d = t_ThermalMacros3d<PetscScalar>;

#endif
