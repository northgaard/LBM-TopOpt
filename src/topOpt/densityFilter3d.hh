#ifndef DENSITYFILTER3D
#define DENSITYFILTER3D

#include "petsc.h"
#include "topOpt/filter.hh"
#include "core/geometry3d.hh"
#include "latticeBoltzmann/LBSolver3d.hh"

class DensityFilter3d : public Filter {

public:
  DensityFilter3d(PetscInt,Box3d,const LBSolver3d&, const AdjointLBSolver&);
  ~DensityFilter3d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
private:
  PetscErrorCode setUp(PetscInt,const Box3d&,DM);
  Mat filterMat;
  Vec denominatorVec;
};

#endif
