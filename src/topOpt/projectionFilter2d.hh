#ifndef PROJECTIONFILTER2D
#define PROJECTIONFILTER2D

#include "petsc.h"
#include "topOpt/filter.hh"
#include "topOpt/densityFilter2d.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

class ProjectionFilter2d : public Filter {

public:
  ProjectionFilter2d(PetscInt,Box2d,const NewLBSolver2d&,const AdjointLBSolver&);
  ~ProjectionFilter2d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
  PetscErrorCode outputIntermediateFields(PetscInt) override;
  void setEta(PetscScalar _e){ eta = _e; }
  void setBeta(PetscScalar _b){ beta = _b; }
  void increaseBeta(PetscScalar inc){ beta *= inc; }
private:
  DensityFilter2d densFilter;
  DM designGrid;
  Vec densityFilteredField;
  Box2d localDomain;
  PetscScalar beta;
  PetscScalar eta;
};

#endif
