#ifndef FINITEDIFFERENCE2D
#define FINITEDIFFERENCE2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"
#include "functionals/LBMacroFunctional.hh"
#include "topOpt/filter.hh"
#include <functional>

class FiniteDifferenceCheck2d {

public:
  using Objective = std::function<PetscErrorCode(NewLBSolver2d&,
                                                 const LBMacroFunctional&,
                                                 PetscScalar*)>;
  using Sensitivities = std::function<PetscErrorCode(NewLBSolver2d&,AdjointLBSolver&,
                                                     const LBMacroFunctional&,
                                                     PetscScalar*,Vec*)>;
  FiniteDifferenceCheck2d(Objective _obj, Sensitivities _sens, PetscScalar _d)
    : computeObjective(_obj), computeSensitivities(_sens), delta(_d){}
  PetscErrorCode checkDomain(NewLBSolver2d&,AdjointLBSolver&,const LBMacroFunctional&,
                             Box2d,Filter* = nullptr,Vec = nullptr);
private:
  Objective computeObjective;
  Sensitivities computeSensitivities;
  PetscScalar delta;
};

#endif

