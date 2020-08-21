#ifndef DENSITYFILTER2D
#define DENSITYFILTER2D

#include "petsc.h"
#include "topOpt/filter.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

class DensityFilter2d : public Filter {

public:
  DensityFilter2d(PetscInt,Box2d,const NewLBSolver2d&,const AdjointLBSolver&);
  ~DensityFilter2d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
protected:
  PetscErrorCode setUp(PetscInt,const Box2d&,DM);
  Mat filterMat;
  Vec denominatorVec;
};

#endif
