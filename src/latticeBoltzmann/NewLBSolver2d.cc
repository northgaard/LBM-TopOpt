#include "NewLBSolver2d.hh"

PetscErrorCode NewLBSolver2d::createPetscObjects(const NewLBSolverInfo2d& info,
                                                 const PetscInt numDOF,
                                                 const PetscInt numMacros,
                                                 const PetscInt numFields)
{
  PetscErrorCode ierr;
  PetscInt stencilWidth = 1;
  DMDAStencilType stencilType = DMDA_STENCIL_BOX;

  ierr = DMDACreate2d(communicator,info.boundaryX,info.boundaryY,stencilType,
                      info.nx,info.ny,PETSC_DECIDE,PETSC_DECIDE,numDOF,stencilWidth,
                      0,0,&latticeGrid);
  CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(latticeGrid,&distributionsGlobal); CHKERRQ(ierr);

  ierr = DMDAGetReducedDMDA(latticeGrid,numMacros,&macroGrid); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(macroGrid,&macrosGlobal); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(macroGrid,&initMacrosGlobal); CHKERRQ(ierr);

  if (numFields){
    ierr = DMDAGetReducedDMDA(latticeGrid,numFields,&materialGrid); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(materialGrid,&materialGlobal); CHKERRQ(ierr);
    ierr = VecSet(materialGlobal,1.); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode NewLBSolver2d::createBoundingBoxes(const NewLBSolverInfo2d& info)
{
  PetscErrorCode ierr;
  PetscInt xs,xm,ys,ym;
  globalBoundingBox = Box2d(0,info.nx-1,0,info.ny-1);
  ierr = DMDAGetCorners(latticeGrid,&xs,&ys,0,&xm,&ym,0); CHKERRQ(ierr);
  localBoundingBox = Box2d(xs,xs+xm-1,ys,ys+ym-1);
  ierr = DMDAGetGhostCorners(latticeGrid,&xs,&ys,0,&xm,&ym,0); CHKERRQ(ierr);
  localBoundingBoxGhosted = Box2d(xs,xs+xm-1,ys,ys+ym-1);

  return 0;
}

PetscErrorCode NewLBSolver2d::setUniformFieldValues(Box2d globalArea, PetscScalar* values)
{
  // TODO: add error check
  PetscErrorCode ierr;
  PetscInt numFields;
  Box2d localArea;
  PetscScalar*** fieldArray;
  PetscFunctionBeginUser;
  boxIntersection(globalArea,localBoundingBox,localArea);
  ierr = DMDAGetInfo(materialGrid,0,0,0,0,0,0,0,&numFields,0,0,0,0,0); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&fieldArray); CHKERRQ(ierr);

  for (auto jj : localArea.yRange){
    for (auto ii : localArea.xRange){
      for (PetscInt ff = 0; ff < numFields; ++ff){
        fieldArray[jj][ii][ff] = values[ff];
      }
    }
  }

  ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&fieldArray); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NewLBSolver2d::addBoundaryConditionRaw(const Box2d boundaryLocation,
                                                      const NewBoundaryDescriptor2d& bcDesc,
                                                      BoundaryOrientation2d orientation)
{
  Box2d localBox;
  PetscErrorCode ierr;
  if (boxIntersection(boundaryLocation,localBoundingBox,localBox)){
    LBBoundaryLoop* loop;
    switch (orientation){
    case North:
      ierr = bcDesc.getNorthBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case South:
      ierr = bcDesc.getSouthBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case West:
      ierr = bcDesc.getWestBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case East:
      ierr = bcDesc.getEastBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case NorthWest:
      ierr = bcDesc.getNorthWestBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case NorthEast:
      ierr = bcDesc.getNorthEastBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case SouthWest:
      ierr = bcDesc.getSouthWestBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    case SouthEast:
      ierr = bcDesc.getSouthEastBoundary(localBox,&loop); CHKERRQ(ierr);
      break;
    }
    boundaryContainer.push_back(loop);
  }
  return 0;
}

PetscErrorCode NewLBSolver2d::addBoundaryCondition(const Box2d boundaryLocation,
                                                   const NewBoundaryDescriptor2d& bcDesc)
{

  PetscErrorCode ierr;
  Box2d tempBox;
  Box2d bcBox;

  if (boxIntersection(boundaryLocation,localBoundingBox,tempBox)){

    Box2dBoundaries boxB = getBoxBoundaries(globalBoundingBox);
    LBBoundaryLoop* loop;

    if (boxIntersection(tempBox,boxB.north,bcBox)){
      ierr = bcDesc.getNorthBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.south,bcBox)){
      ierr = bcDesc.getSouthBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.west,bcBox)){
      ierr = bcDesc.getWestBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.east,bcBox)){
      ierr = bcDesc.getEastBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northWest,bcBox)){
      ierr = bcDesc.getNorthWestBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northEast,bcBox)){
      ierr = bcDesc.getNorthEastBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southWest,bcBox)){
      ierr = bcDesc.getSouthWestBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southEast,bcBox)){
      ierr = bcDesc.getSouthEastBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
  }

  return 0;

}
