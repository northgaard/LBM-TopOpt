#include "obstacleThermalLBSolver2d.hh"

PetscErrorCode ObstacleThermalLBSolver2d::collide()
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

PetscErrorCode ObstacleThermalLBSolver2d::streamAndCollide()
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

PetscErrorCode ObstacleThermalLBSolver2d::
uniformMacroInitialization(const ThermalMacros2d& values)
{

  PetscErrorCode ierr;
  ThermalMacros2d** initMac;

  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,&initMac); CHKERRQ(ierr);

  for (auto jj : localBoundingBox.yRange){
    for (auto ii : localBoundingBox.xRange){
      initMac[jj][ii].rho = values.rho;
      initMac[jj][ii].ux = values.ux;
      initMac[jj][ii].uy = values.uy;
      initMac[jj][ii].T = values.T;
    }
  }

  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,&initMac); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ObstacleThermalLBSolver2d::
getInitialMacroArray(const Box2d& inputBox,
		     Box2d& outLocalBox,
		     ThermalMacros2d*** theArray)
{

  PetscErrorCode ierr;
  boxIntersection(inputBox,localBoundingBox,outLocalBox);
  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,theArray); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ObstacleThermalLBSolver2d::
restoreInitialMacroArray(const Box2d& inputBox,
			 Box2d& outLocalBox,
			 ThermalMacros2d*** theArray)
{

  PetscErrorCode ierr;
  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,theArray); CHKERRQ(ierr);

  return 0;

}

