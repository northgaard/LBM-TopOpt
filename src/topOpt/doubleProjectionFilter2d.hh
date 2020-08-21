#ifndef DOUBLEPROJECTIONFILER2D
#define DOUBLEPROJECTIONFILER2D

#include "petsc.h"
#include "topOpt/filter.hh"
#include "topOpt/projectionFilter2d.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

class DoubleProjectionFilter2d : public Filter {

public:
  DoubleProjectionFilter2d(PetscInt,PetscInt,Box2d,
                           const NewLBSolver2d&, const AdjointLBSolver&);
  ~DoubleProjectionFilter2d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
  PetscErrorCode outputIntermediateFields(PetscInt) override;
  void setFirstEta(PetscScalar _e){ firstProjection.setEta(_e); }
  void setSecondEta(PetscScalar _e){ secondProjection.setEta(_e); }
  void setFirstBeta(PetscScalar _b){ firstProjection.setBeta(_b); }
  void setSecondBeta(PetscScalar _b){ secondProjection.setBeta(_b); }
  void increaseBeta(PetscScalar inc)
  {
    firstProjection.increaseBeta(inc);
    secondProjection.increaseBeta(inc);
  }
private:
  ProjectionFilter2d firstProjection;
  ProjectionFilter2d secondProjection;
  Vec firstProjectedField;
};

#endif
