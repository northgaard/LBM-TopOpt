#include "LBSolverBase2d.hh"

PetscErrorCode LBSolverBase2d::stream()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** fcol;

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  streamFunc(fdist,fcol);
  ++currentTimestep;
  boundaryComputations(currentTimestep,fdist);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsLocal,&fcol); CHKERRQ(ierr);

  return 0;
  
}

void LBSolverBase2d::boundaryComputations(PetscInt timestep, PetscScalar*** fdist)
{
  for (const auto& fun : boundaryContainer){
    fun(timestep,fdist);
  }
}

PetscErrorCode LBSolverBase2d::initializeDistributions()
{

  PetscErrorCode ierr;
  PetscScalar*** fdist;
  PetscScalar*** mac;

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,initMacrosGlobal,&mac); CHKERRQ(ierr);

  initialization(fdist,mac);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,initMacrosGlobal,&mac); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
			      distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
			    distributionsLocal); CHKERRQ(ierr);
  ierr = VecCopy(initMacrosGlobal,macrosGlobal); CHKERRQ(ierr);

  currentTimestep = 0;

  return 0;

}

PetscErrorCode LBSolverBase2d::createLatticeObjects(PetscInt nx, PetscInt ny,
						    PetscInt numDOF,
						    PetscDMBoundary boundaryX,
						    PetscDMBoundary boundaryY)
{

  PetscErrorCode ierr;
  PetscInt stencilWidth = 1;
  DMDAStencilType stencilType = DMDA_STENCIL_BOX;

  // Create lattice DMDA and vectors
  ierr = DMDACreate2d(PETSC_COMM_WORLD,boundaryX,boundaryY,stencilType,
		      nx,ny,PETSC_DECIDE,PETSC_DECIDE,numDOF,stencilWidth,
		      0,0,&latticeGrid);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) latticeGrid);
  CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(latticeGrid,&distributionsGlobal);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) distributionsGlobal);
  CHKERRQ(ierr);
  ierr = DMCreateLocalVector(latticeGrid,&distributionsLocal);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) distributionsLocal);
  CHKERRQ(ierr);

  return 0;
  
}

PetscErrorCode LBSolverBase2d::createMacroObjects(PetscInt numMac)
{

  PetscErrorCode ierr;

  ierr = DMDAGetReducedDMDA(latticeGrid,numMac,&macroGrid);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) macroGrid);
  CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(macroGrid,&initMacrosGlobal);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) initMacrosGlobal);
  CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(macroGrid,&macrosGlobal);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) macrosGlobal);
  CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBSolverBase2d::createBoundingBoxes(PetscInt nx, PetscInt ny)
{

  PetscErrorCode ierr;
  PetscInt xs,xm,ys,ym;

  globalBoundingBox = Box2d(0,nx-1,0,ny-1);
  ierr = DMDAGetCorners(latticeGrid,&xs,&ys,0,&xm,&ym,0);
  CHKERRQ(ierr);
  localBoundingBox = Box2d(xs,xs+xm-1,ys,ys+ym-1);

  ierr = DMDAGetGhostCorners(latticeGrid,&xs,&ys,0,&xm,&ym,0);
  CHKERRQ(ierr);
  localBoundingBoxGhosted = Box2d(xs,xs+xm-1,ys,ys+ym-1);

  return 0;

}

PetscErrorCode LBSolverBase2d::outputMacros()
{

  PetscErrorCode ierr;

  PetscViewer view;
  char output[50];
  sprintf(output,"Macro_data_iter_%i.vts",currentTimestep);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,output,FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(macrosGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBSolverBase2d::outputDistributions()
{

  PetscErrorCode ierr;

  PetscViewer view;
  char output[50];
  sprintf(output,"Distribution_data_iter_%i.vts",currentTimestep);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,output,FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(distributionsGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBSolverBase2d::
computeMacroObjective(const std::unique_ptr<ObjectiveFunction2d>& theFunc,
		      PetscScalar& value)
{

  PetscErrorCode ierr;
  PetscScalar*** mac;

  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  theFunc->evaluate(mac,value);

  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBSolverBase2d::
addBoundaryCondition(const Box2d boundaryLocation,
		     const std::unique_ptr<BoundaryDescriptor2d>& bcDesc)
{

  PetscErrorCode ierr;
  Box2d tempBox;
  Box2d bcBox;

  if (boxIntersection(boundaryLocation,localBoundingBox,tempBox)){

    Box2dBoundaries boxB = getBoxBoundaries(globalBoundingBox);

    if (boxIntersection(tempBox,boxB.north,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getNorthBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.south,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getSouthBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.west,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getWestBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.east,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getEastBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.northWest,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getNorthWestBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.northEast,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getNorthEastBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.southWest,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getSouthWestBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
    if (boxIntersection(tempBox,boxB.southEast,bcBox)){
      BoundaryFunction func;
      ierr = bcDesc->getSouthEastBoundary(bcBox,func); CHKERRQ(ierr);
      boundaryContainer.push_back(func);
    }
  }

  return 0;

}
