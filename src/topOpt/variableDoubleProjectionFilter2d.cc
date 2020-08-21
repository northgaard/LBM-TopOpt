#include "variableDoubleProjectionFilter2d.hh"

VariableDoubleProjectionFilter2d::
VariableDoubleProjectionFilter2d(PetscInt firstFilterRadius,
                                 PetscInt secondFilterRadius,
                                 Box2d filterDomain,
                                 const NewLBSolver2d& solver,
                                 const AdjointLBSolver& adjSolver)
  : firstProjection(firstFilterRadius,filterDomain,solver,adjSolver),
    secondProjection(secondFilterRadius,filterDomain,solver,adjSolver),
    firstProjectedField(nullptr)
{
  PetscErrorCode ierr;
  ierr = DMCreateGlobalVector(adjSolver.sensitivityGrid,&firstProjectedField);
  CHKERRABORT(solver.communicator,ierr);
}

VariableDoubleProjectionFilter2d::~VariableDoubleProjectionFilter2d()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  if (firstProjectedField){
    PetscObjectGetComm((PetscObject) firstProjectedField, &comm);
    ierr = VecDestroy(&firstProjectedField);
    CHKERRABORT(comm,ierr);
  }
}

PetscErrorCode VariableDoubleProjectionFilter2d::
setSecondEtaValuesFromVector(Vec theVec)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = secondProjection.setEtaValuesFromVector(theVec); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableDoubleProjectionFilter2d::
filterDesign(Vec designField, Vec projectedField)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = firstProjection.filterDesign(designField,firstProjectedField); CHKERRQ(ierr);
  ierr = secondProjection.filterDesign(firstProjectedField,projectedField); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableDoubleProjectionFilter2d::
filterSensitivities(Vec sensitivityField, Vec* constraintSensitivityFields,
                    PetscInt numConstraints)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = secondProjection.filterSensitivities(sensitivityField,constraintSensitivityFields,
                                              numConstraints); CHKERRQ(ierr);
  ierr = firstProjection.filterSensitivities(sensitivityField,constraintSensitivityFields,
                                             numConstraints); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableDoubleProjectionFilter2d::
outputSecondEtaValues(const std::string& outputFolder, const std::string& output)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = secondProjection.outputEtaValues(outputFolder,output); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VariableDoubleProjectionFilter2d::
meanSecondEtaValue(PetscScalar* value)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = secondProjection.meanEtaValue(value); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// TODO: implement this, might have to change the interface a bit
PetscErrorCode VariableDoubleProjectionFilter2d::outputIntermediateFields(PetscInt)
{
  return 0;
}
