#include "LBSolver2d.hh"

PetscErrorCode LBSolver2d::initializeSolver(const LBSolverInfo2d& info)
{

  PetscErrorCode ierr;
  ierr = createLatticeObjects(info.nx,info.ny,info.latticeDOF,
			      info.boundaryX,info.boundaryY);
  CHKERRQ(ierr);
  ierr = createMacroObjects(info.macroDOF); CHKERRQ(ierr);
  ierr = createBoundingBoxes(info.nx,info.ny); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBSolver2d::collide()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** mac;

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  collisionFunc(fdist,mac);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
			      distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
			    distributionsLocal); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBSolver2d::streamAndCollide()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** fcol;
  PetscScalar*** mac;

  /* Stream */

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  streamFunc(fdist,fcol);
  ++currentTimestep;
  boundaryComputations(currentTimestep,fdist);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  /* Collide */

  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  collisionFunc(fdist,mac);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
			      distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
			    distributionsLocal); CHKERRQ(ierr);

  return 0;
  
}

/* Should test if the Latt "single loop" is faster for this! */

PetscErrorCode LBSolver2d::collideAndStream()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** fcol;
  PetscScalar*** mac;

  /* Collide */
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  collisionFunc(fdist,mac);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,distributionsGlobal,&mac); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
                              distributionsLocal); CHKERRQ(ierr); 
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
                            distributionsLocal); CHKERRQ(ierr); 

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  streamFunc(fdist,fcol);
  ++currentTimestep;
  boundaryComputations(currentTimestep,fdist);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  return 0;

}
