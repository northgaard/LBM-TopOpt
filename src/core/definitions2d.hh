#ifndef DEFINITIONS2D
#define DEFINITIONS2D

#include "petsc.h"

template <class Real>
struct t_IsothermalMacros2d {
  Real rho;
  Real ux;
  Real uy;
};

using IsothermalMacros2d = t_IsothermalMacros2d<PetscScalar>;

template <class Real>
struct t_ThermalMacros2d {
  Real rho;
  Real ux;
  Real uy;
  Real T;
};

using ThermalMacros2d = t_ThermalMacros2d<PetscScalar>;

enum BoundaryOrientation2d {North,South,West,East,NorthWest,NorthEast,SouthWest,SouthEast};

#endif
