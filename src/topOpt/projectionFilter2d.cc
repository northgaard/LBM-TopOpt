#include "projectionFilter2d.hh"
#include <cmath>

ProjectionFilter2d::ProjectionFilter2d(PetscInt filterRadius, Box2d filterDomain,
                                       const NewLBSolver2d& solver,
                                       const AdjointLBSolver& adjSolver)
  : densFilter(filterRadius,filterDomain,solver,adjSolver), beta(1.), eta(0.5)
{
  PetscErrorCode ierr;
  boxIntersection(filterDomain,solver.getLocalBoundingBox(),localDomain);
  designGrid = adjSolver.sensitivityGrid;
  ierr = DMCreateGlobalVector(designGrid,&densityFilteredField);
  CHKERRABORT(solver.communicator,ierr);
}

ProjectionFilter2d::~ProjectionFilter2d()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  if (densityFilteredField){
    PetscObjectGetComm((PetscObject) densityFilteredField, &comm);
    ierr = VecDestroy(&densityFilteredField);
    CHKERRABORT(comm,ierr);
  }
}

PetscErrorCode ProjectionFilter2d::filterDesign(Vec designField, Vec projectedField)
{
  PetscErrorCode ierr;
  PetscScalar **filteredArr, **projectedArr;
  PetscFunctionBeginUser;
  ierr = densFilter.filterDesign(designField,densityFilteredField); CHKERRQ(ierr);
  ierr = VecCopy(densityFilteredField,projectedField); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,densityFilteredField,&filteredArr); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,projectedField,&projectedArr); CHKERRQ(ierr);

  for (auto jj : localDomain.yRange){
    for (auto ii : localDomain.xRange){
      projectedArr[jj][ii] = (tanh(beta*eta) + tanh(beta*(filteredArr[jj][ii] - eta)))
        / (tanh(beta*eta) + tanh(beta*(1. - eta)));
    }
  }

  ierr = DMDAVecRestoreArray(designGrid,densityFilteredField,&filteredArr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(designGrid,projectedField,&projectedArr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ProjectionFilter2d::filterSensitivities(Vec sensitivityField,
                                                       Vec* constraintSensitivityFields,
                                                       PetscInt numConstraints)
{
  PetscErrorCode ierr;
  PetscScalar** sens, **densityFilterArr;
  PetscScalar*** constraintSens;
  PetscScalar projDiff;
  PetscFunctionBeginUser;
  ierr = PetscMalloc1(numConstraints,&constraintSens); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,densityFilteredField,&densityFilterArr); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(designGrid,sensitivityField,&sens); CHKERRQ(ierr);
  for (PetscInt cc = 0; cc < numConstraints; ++cc){
    ierr = DMDAVecGetArray(designGrid,constraintSensitivityFields[cc],
                           &(constraintSens[cc])); CHKERRQ(ierr);
  }

  for (auto jj : localDomain.yRange){
    for (auto ii : localDomain.xRange){
      projDiff = ((1.0 - pow(tanh(beta*(densityFilterArr[jj][ii] - eta)),2.0))*beta)
        / (tanh(beta*eta) + tanh(beta*(1.0 - eta)));
      sens[jj][ii] *= projDiff;
      for (PetscInt cc = 0; cc < numConstraints; ++cc){
        constraintSens[cc][jj][ii] *= projDiff;
      }
    }
  }

  ierr = DMDAVecRestoreArray(designGrid,densityFilteredField,&densityFilterArr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(designGrid,sensitivityField,&sens); CHKERRQ(ierr);
  for (PetscInt cc = 0; cc < numConstraints; ++cc){
    ierr = DMDAVecRestoreArray(designGrid,constraintSensitivityFields[cc],
                               &(constraintSens[cc])); CHKERRQ(ierr);
  }
  ierr = densFilter.filterSensitivities(sensitivityField,constraintSensitivityFields,
                                        numConstraints); CHKERRQ(ierr);
  ierr = PetscFree(constraintSens); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ProjectionFilter2d::outputIntermediateFields(PetscInt numIteration)
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
