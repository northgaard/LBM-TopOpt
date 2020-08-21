#include "variableProjectionFilter2d.hh"
#include <cmath>

VariableProjectionFilter2d::VariableProjectionFilter2d(PetscInt filterRadius,
                                                       Box2d filterDomain,
                                                       const NewLBSolver2d& solver,
                                                       const AdjointLBSolver& adjSolver)
  : densFilter(filterRadius,filterDomain,solver,adjSolver), beta(1.), etaValues(nullptr)
{
  PetscErrorCode ierr;
  boxIntersection(filterDomain,solver.getLocalBoundingBox(),localDomain);
  designGrid = adjSolver.sensitivityGrid;
  ierr = DMCreateGlobalVector(designGrid,&densityFilteredField);
  CHKERRABORT(solver.communicator,ierr);
  ierr = DMCreateGlobalVector(designGrid,&etaValues);
  CHKERRABORT(solver.communicator,ierr);
  ierr = VecSet(etaValues,0.5);
  CHKERRABORT(solver.communicator,ierr);
  globalDomainSize = filterDomain.getNx()*filterDomain.getNy();
}

VariableProjectionFilter2d::~VariableProjectionFilter2d()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  if (densityFilteredField){
    PetscObjectGetComm((PetscObject) densityFilteredField, &comm);
    ierr = VecDestroy(&densityFilteredField);
    CHKERRABORT(comm,ierr);
  }
  if (etaValues){
    PetscObjectGetComm((PetscObject) etaValues, &comm);
    ierr = VecDestroy(&etaValues);
    CHKERRABORT(comm,ierr);
  }
}


PetscErrorCode VariableProjectionFilter2d::filterDesign(Vec designField,
                                                        Vec projectedField)
{
  PetscErrorCode ierr;
  PetscScalar **filteredArr, **projectedArr, **eta;
  PetscFunctionBeginUser;
  ierr = densFilter.filterDesign(designField,densityFilteredField); CHKERRQ(ierr);
  ierr = VecCopy(densityFilteredField,projectedField); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,densityFilteredField,&filteredArr); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,projectedField,&projectedArr); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,etaValues,&eta); CHKERRQ(ierr);

  for (auto jj : localDomain.yRange){
    for (auto ii : localDomain.xRange){
      projectedArr[jj][ii] =
        (tanh(beta*eta[jj][ii]) + tanh(beta*(filteredArr[jj][ii] - eta[jj][ii])))
        / (tanh(beta*eta[jj][ii]) + tanh(beta*(1. - eta[jj][ii])));
    }
  }

  ierr = DMDAVecRestoreArray(designGrid,densityFilteredField,&filteredArr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(designGrid,projectedField,&projectedArr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(designGrid,etaValues,&eta); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableProjectionFilter2d::filterSensitivities(Vec sensitivityField,
                                                               Vec* constraintSensitivityFields,
                                                               PetscInt numConstraints)
{
  PetscErrorCode ierr;
  PetscScalar** sens, **densityFilterArr, **eta;
  PetscScalar*** constraintSens;
  PetscScalar projDiff;
  PetscFunctionBeginUser;
  ierr = PetscMalloc1(numConstraints,&constraintSens); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,densityFilteredField,&densityFilterArr); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,sensitivityField,&sens); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,etaValues,&eta); CHKERRQ(ierr);
  for (PetscInt cc = 0; cc < numConstraints; ++cc){
    ierr = DMDAVecGetArray(designGrid,constraintSensitivityFields[cc],
                           &(constraintSens[cc])); CHKERRQ(ierr);
  }

  for (auto jj : localDomain.yRange){
    for (auto ii : localDomain.xRange){
      projDiff = ((1.0 - pow(tanh(beta*(densityFilterArr[jj][ii] - eta[jj][ii])),2.0))*beta)
        / (tanh(beta*eta[jj][ii]) + tanh(beta*(1.0 - eta[jj][ii])));
      sens[jj][ii] *= projDiff;
      for (PetscInt cc = 0; cc < numConstraints; ++cc){
        constraintSens[cc][jj][ii] *= projDiff;
      }
    }
  }

  ierr = DMDAVecRestoreArray(designGrid,densityFilteredField,&densityFilterArr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(designGrid,sensitivityField,&sens); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(designGrid,etaValues,&eta); CHKERRQ(ierr);
  for (PetscInt cc = 0; cc < numConstraints; ++cc){
    ierr = DMDAVecRestoreArray(designGrid,constraintSensitivityFields[cc],
                               &(constraintSens[cc])); CHKERRQ(ierr);
  }
  ierr = densFilter.filterSensitivities(sensitivityField,constraintSensitivityFields,
                                        numConstraints); CHKERRQ(ierr);
  ierr = PetscFree(constraintSens); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableProjectionFilter2d::outputIntermediateFields(PetscInt numIteration)
{
  PetscErrorCode ierr;
  PetscViewer view;
  char output[50];
  MPI_Comm comm;
  PetscFunctionBeginUser;
  sprintf(output,"Density_filtered_field_iter_%i.vts",numIteration);
  PetscObjectGetComm((PetscObject) densityFilteredField, &comm);
  ierr = PetscViewerVTKOpen(comm,output,FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(densityFilteredField,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableProjectionFilter2d::setEtaValuesFromVector(Vec theVec)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecCopy(theVec,etaValues); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableProjectionFilter2d::outputEtaValues(const std::string& folder,
                                                           const std::string& name)
{
  PetscErrorCode ierr;
  PetscViewer view;
  MPI_Comm communicator;
  PetscFunctionBeginUser;
  std::string output;
  output = folder + name;
  ierr = PetscObjectGetComm((PetscObject) etaValues,&communicator); CHKERRQ(ierr);
  ierr = PetscViewerVTKOpen(communicator,output.c_str(),FILE_MODE_WRITE,&view); CHKERRQ(ierr);
  ierr = VecView(etaValues,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableProjectionFilter2d::meanEtaValue(PetscScalar* value)
{
  PetscErrorCode ierr;
  PetscScalar **etaArr;
  PetscScalar lsum = 0.;
  MPI_Comm communicator;
  PetscFunctionBeginUser;
  ierr = PetscObjectGetComm((PetscObject) etaValues,&communicator); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,etaValues,&etaArr); CHKERRQ(ierr);

  for (auto jj : localDomain.yRange){
    for (auto ii : localDomain.xRange){
      lsum += etaArr[jj][ii];
    }
  }
  MPI_Allreduce(&lsum,value,1,MPIU_SCALAR,MPI_SUM,communicator);
  *value /= globalDomainSize;

  ierr = DMDAVecRestoreArray(designGrid,etaValues,&etaArr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
