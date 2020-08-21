#ifndef VARIABLEPROJECTIONFILTER2D
#define VARIABLEPROJECTIONFILTER2D

#include "petsc.h"
#include "topOpt/filter.hh"
#include "topOpt/densityFilter2d.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"
#include <string>

class VariableProjectionFilter2d : public Filter {

public:
  VariableProjectionFilter2d(PetscInt,Box2d,const NewLBSolver2d&,const AdjointLBSolver&);
  ~VariableProjectionFilter2d();
  PetscErrorCode filterDesign(Vec,Vec) override;
  PetscErrorCode filterSensitivities(Vec,Vec*,PetscInt) override;
  PetscErrorCode outputIntermediateFields(PetscInt) override;
  void setBeta(PetscScalar _b){ beta = _b; }
  void increaseBeta(PetscScalar inc){ beta *= inc; }
  PetscErrorCode setEtaValuesFromVector(Vec);
  template <typename Callable>
  PetscErrorCode setEtaValuesFromFunction(const Callable& call)
  {
    PetscErrorCode ierr;
    PetscScalar **etaArr;
    PetscFunctionBeginUser;
    ierr = DMDAVecGetArray(designGrid,etaValues,&etaArr); CHKERRQ(ierr);

    for (auto jj : localDomain.yRange){
      for (auto ii : localDomain.xRange){
        etaArr[jj][ii] = call(ii,jj);
      }
    }

    ierr = DMDAVecRestoreArray(designGrid,etaValues,&etaArr); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode outputEtaValues(const std::string&,const std::string&);
  PetscErrorCode meanEtaValue(PetscScalar*);
private:
  DensityFilter2d densFilter;
  DM designGrid;
  Vec densityFilteredField;
  Box2d localDomain;
  PetscScalar beta;
  Vec etaValues;
  PetscInt globalDomainSize;
};

#endif
