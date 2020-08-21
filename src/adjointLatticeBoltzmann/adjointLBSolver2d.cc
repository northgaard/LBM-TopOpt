#include "adjointLBSolver2d.hh"

PetscErrorCode AdjointLBSolver2d::adjointCollide(Vec distributionsGlobal,
						 Vec obstacleGlobal,
						 Vec sensitivity,
						 ObjectiveFunction2d* obj)
{

  PetscErrorCode ierr;
  PetscScalar*** adj;
  PetscScalar*** fdist;
  PetscScalar** obst;
  PetscScalar** os;

  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(obstacleGrid,sensitivity,&os); CHKERRQ(ierr);

  adjointCollisionFunc(adj,fdist,obst,os);
  obj->adjointStateSource(curTimestep,adj,fdist);
  boundaryComputations(curTimestep,adj,fdist);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(obstacleGrid,obstacleGlobal,&obst); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(obstacleGrid,sensitivity,&os); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,adjointGlobal,INSERT_VALUES,
			      adjointLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,adjointGlobal,INSERT_VALUES,
			    adjointLocal); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode AdjointLBSolver2d::adjointStream()
{

  PetscErrorCode ierr;
  PetscScalar*** adj;
  PetscScalar*** adjPrev;

  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointLocal,&adjPrev); CHKERRQ(ierr);

  --curTimestep;
  adjointStreamFunc(adj,adjPrev);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointLocal,&adjPrev); CHKERRQ(ierr);

  return 0;

}

void AdjointLBSolver2d::boundaryComputations(PetscInt timestep,
					     PetscScalar*** adj,
					     PetscScalar*** fdist)
{
  for (const auto& fun : adjointBoundaryContainer){
    fun(timestep,adj,fdist);
  }
}

PetscErrorCode AdjointLBSolver2d::
initializeFromForwardSolver(const ObstacleLBSolver2d& solver)
{

  latticeGrid = solver.latticeGrid;
  macroGrid = solver.macroGrid;
  obstacleGrid = solver.obstacleGrid;

  globalBoundingBox = solver.getBoundingBox();
  localBoundingBox = solver.getLocalBoundingBox();
  localBoundingBoxGhosted = solver.getLocalBoundingBoxGhosted();

  PetscErrorCode ierr;

  ierr = createAdjointVectors(); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode AdjointLBSolver2d::createAdjointVectors()
{

  PetscErrorCode ierr;
  ierr = DMCreateGlobalVector(latticeGrid,&adjointGlobal); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) adjointGlobal); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(latticeGrid,&adjointLocal); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) adjointLocal); CHKERRQ(ierr);

  ierr = VecSet(adjointGlobal,0.); CHKERRQ(ierr);
  ierr = VecSet(adjointLocal,0.); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode AdjointLBSolver2d::
addBoundaryCondition(const Box2d boundaryLocation,
		     const std::unique_ptr<BoundaryDescriptor2d>& bcDesc)
{

  PetscErrorCode ierr;
  Box2d tempBox;
  Box2d bcBox;

  if (boxIntersection(boundaryLocation,localBoundingBox,tempBox)){

    Box2dBoundaries boxB = getBoxBoundaries(globalBoundingBox);

    if (boxIntersection(tempBox,boxB.north,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointNorthBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.south,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointSouthBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.west,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointWestBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.east,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointEastBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.northWest,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointNorthWestBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.northEast,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointNorthEastBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.southWest,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointSouthWestBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.southEast,bcBox)){
      AdjointBoundaryFunction func;
      ierr = bcDesc->getAdjointSouthEastBoundary(bcBox,func); CHKERRQ(ierr);
      adjointBoundaryContainer.push_back(func);
    }
  }

  return 0;

}

PetscErrorCode
addBoundaryCondition(const Box2d boundaryLocation,
		     const std::unique_ptr<BoundaryDescriptor2d>& bcDesc,
		     ObstacleLBSolver2d& forward,
		     AdjointLBSolver2d& adjoint)
{

  PetscErrorCode ierr;
  ierr = forward.addBoundaryCondition(boundaryLocation,bcDesc); CHKERRQ(ierr);
  ierr = adjoint.addBoundaryCondition(boundaryLocation,bcDesc); CHKERRQ(ierr);

  return 0;

}
