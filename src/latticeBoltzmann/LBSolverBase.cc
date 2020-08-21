#include "LBSolverBase.hh"

LBSolverBase::LBSolverBase() : distributionsGlobal(0), distributionsLocal(0),
                               macrosGlobal(0), initMacrosGlobal(0){}
LBSolverBase::~LBSolverBase(){} 

PetscErrorCode LBSolverBase::setDistributions(Vec dist)
{
  // PetscErrorCode ierr;
  // ierr = VecCopy(dist,distributionsGlobal); CHKERRQ(ierr);

  return 0;
}
