#ifndef VARIABLEDOUBLEPROJECTIONFILTER2D
#define VARIABLEDOUBLEPROJECTIONFILTER2D

#include "petsc.h"
#include "topOpt/filter.hh"
#include "topOpt/projectionFilter2d.hh"
#include "topOpt/variableProjectionFilter2d.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

class VariableDoubleProjectionFilter2d : public Filter {

public:
  VariableDoubleProjectionFilter2d(PetscInt,PetscInt,Box2d,
                                   const NewLBSolver2d&, const AdjointLBSolver&);
  ~VariableDoubleProjectionFilter2d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
  PetscErrorCode outputIntermediateFields(PetscInt) override;
  void setFirstEta(PetscScalar _e){ firstProjection.setEta(_e); }
  void setFirstBeta(PetscScalar _b){ firstProjection.setBeta(_b); }
  void setSecondBeta(PetscScalar _b){ secondProjection.setBeta(_b); }
  void increaseBeta(PetscScalar inc)
  {
    firstProjection.increaseBeta(inc);
    secondProjection.increaseBeta(inc);
  }
  PetscErrorCode setSecondEtaValuesFromVector(Vec);
  template <typename Callable>
  PetscErrorCode setSecondEtaValuesFromFunction(const Callable& call)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = secondProjection.setEtaValuesFromFunction(call); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode outputSecondEtaValues(const std::string&,
                                       const std::string&);
  PetscErrorCode meanSecondEtaValue(PetscScalar*);
private:
  ProjectionFilter2d firstProjection;
  VariableProjectionFilter2d secondProjection;
  Vec firstProjectedField;
};

#endif
