#include "obstacleLBSolverBase2d.hh"

PetscErrorCode ObstacleLBSolverBase2d::createObstacleObjects(PetscInt numDOF)
{

  PetscErrorCode ierr;

  ierr = DMDAGetReducedDMDA(latticeGrid,numDOF,&obstacleGrid); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) obstacleGrid);
  CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(obstacleGrid,&obstacleGlobal); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) obstacleGlobal);
  CHKERRQ(ierr);
  ierr = VecSet(obstacleGlobal,1.); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ObstacleLBSolverBase2d::addGeometricObstacle(const Box2d& theBox,
							    PetscScalar value)
{

  PetscErrorCode ierr;
  PetscScalar** obst;
  Box2d lb;

  ierr = DMDAVecGetArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  if (boxIntersection(theBox,localBoundingBox,lb)){
    for (auto jj : lb.yRange){
      for (auto ii : lb.xRange){
	obst[jj][ii] = value;
      }
    }
  }

  ierr = DMDAVecRestoreArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode ObstacleLBSolverBase2d::outputDomain(std::string output)
{

  PetscErrorCode ierr;
  PetscViewer view;

  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,output.c_str(),FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(obstacleGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  
  return 0;

}

PetscErrorCode ObstacleLBSolverBase2d::setDomain(Vec newDomain)
{

  PetscErrorCode ierr;
  ierr = VecCopy(newDomain,obstacleGlobal); CHKERRQ(ierr);

  return 0;

}
