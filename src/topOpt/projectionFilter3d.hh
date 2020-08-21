#ifndef PROJECTIONFILTER3D
#define PROJECTIONFILTER3D

#include "petsc.h"
#include "topOpt/densityFilter3d.hh"
#include "core/geometry3d.hh"
#include "latticeBoltzmann/LBSolver3d.hh"

class ProjectionFilter3d : public Filter {

public:
  ProjectionFilter3d(PetscInt,Box3d,const LBSolver3d&,const AdjointLBSolver&);
  ~ProjectionFilter3d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
  PetscErrorCode outputIntermediateFields(PetscInt) override;
  void setEta(PetscScalar _e){ eta = _e; }
  void setBeta(PetscScalar _b){ beta = _b; }
  void increaseBeta(PetscScalar inc){ beta *= inc; }
private:
  DensityFilter3d densFilter;
  DM designGrid;
  Vec densityFilteredField;
  Box3d localDomain;
  PetscScalar beta;
  PetscScalar eta;
};

#endif
