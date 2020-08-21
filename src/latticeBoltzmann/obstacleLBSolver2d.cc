#include "obstacleLBSolver2d.hh"

PetscErrorCode ObstacleLBSolver2d::initializeSolver(const ObstacleLBSolverInfo2d& info)
{

  PetscErrorCode ierr;
  ierr = createLatticeObjects(info.nx,info.ny,info.latticeDOF,
			      info.boundaryX,info.boundaryY);
  CHKERRQ(ierr);
  ierr = createMacroObjects(info.macroDOF); CHKERRQ(ierr);
  ierr = createObstacleObjects(info.obstacleDOF); CHKERRQ(ierr);
  ierr = createBoundingBoxes(info.nx,info.ny); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ObstacleLBSolver2d::collide()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** mac;
  PetscScalar** obst;

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  collisionFunc(fdist,mac,obst);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
			      distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
			    distributionsLocal); CHKERRQ(ierr);

  return 0;
  
}

PetscErrorCode ObstacleLBSolver2d::streamAndCollide()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** fcol;
  PetscScalar*** mac;
  PetscScalar** obst;

  /* Stream */

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  streamFunc(fdist,fcol);
  ++currentTimestep;
  boundaryComputations(currentTimestep,fdist);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  /* Collide */

  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  collisionFunc(fdist,mac,obst);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
			      distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
			    distributionsLocal); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ObstacleLBSolver2d::collideAndStream()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** fcol;
  PetscScalar*** mac;
  PetscScalar** obst;

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  collisionFunc(fdist,mac,obst);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

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
