#include "LBSolver3d.hh"

PetscErrorCode LBSolver3d::createPetscObjects(const LBSolverInfo3d& info,
                                              const PetscInt numDOF,
                                              const PetscInt numMacros,
                                              const PetscInt numFields)
{
  PetscErrorCode ierr;
  PetscInt stencilWidth = 1;
  DMDAStencilType stencilType = DMDA_STENCIL_BOX;
  PetscFunctionBeginUser;
  ierr = DMDACreate3d(PETSC_COMM_WORLD,info.boundaryX,info.boundaryY,info.boundaryZ,
                      stencilType,info.nx,info.ny,info.nz,PETSC_DECIDE,PETSC_DECIDE,
                      PETSC_DECIDE,numDOF,stencilWidth,0,0,0,&latticeGrid);
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
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver3d::createBoundingBoxes(const LBSolverInfo3d& info)
{
  PetscErrorCode ierr;
  PetscInt xs,xm,ys,ym,zs,zm;
  PetscFunctionBeginUser;
  globalBoundingBox = Box3d(0,info.nx-1,0,info.ny-1,0,info.nz-1);
  ierr = DMDAGetCorners(latticeGrid,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  localBoundingBox = Box3d(xs,xs+xm-1,ys,ys+ym-1,zs,zs+zm-1);
  ierr = DMDAGetGhostCorners(latticeGrid,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  localBoundingBoxGhosted = Box3d(xs,xs+xm-1,ys,ys+ym-1,zs,zs+zm-1);
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver3d::setUniformFieldValues(Box3d globalArea, PetscScalar* values)
{
  PetscErrorCode ierr;
  PetscInt numFields;
  Box3d localArea;
  PetscScalar**** fieldArray;
  PetscFunctionBeginUser;
  boxIntersection(globalArea,localBoundingBox,localArea);
  ierr = DMDAGetInfo(materialGrid,0,0,0,0,0,0,0,&numFields,0,0,0,0,0); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&fieldArray); CHKERRQ(ierr);

  for (auto kk : localArea.zRange){
    for (auto jj : localArea.yRange){
      for (auto ii : localArea.xRange){
        for (PetscInt ff = 0; ff < numFields; ++ff){
          fieldArray[kk][jj][ii][ff] = values[ff];
        }
      }
    }
  }

  ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&fieldArray); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver3d::addBoundaryCondition(const Box3d boundaryLocation,
                                                const BoundaryDescriptor3d& bcDesc)
{
  PetscErrorCode ierr;
  Box3d tempBox;
  Box3d bcBox;
  PetscFunctionBeginUser;
  if (boxIntersection(boundaryLocation,localBoundingBox,tempBox)){
    Box3dBoundaries boxB = getBoxBoundaries(globalBoundingBox);
    LBBoundaryLoop* loop;
    /* Check faces */
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
    if (boxIntersection(tempBox,boxB.front,bcBox)){
      ierr = bcDesc.getFrontBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.back,bcBox)){
      ierr = bcDesc.getBackBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    /* Check edges */
    if (boxIntersection(tempBox,boxB.northWestEdge,bcBox)){
      ierr = bcDesc.getNorthWestEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northEastEdge,bcBox)){
      ierr = bcDesc.getNorthEastEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northFrontEdge,bcBox)){
      ierr = bcDesc.getNorthFrontEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northBackEdge,bcBox)){
      ierr = bcDesc.getNorthBackEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southWestEdge,bcBox)){
      ierr = bcDesc.getSouthWestEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southEastEdge,bcBox)){
      ierr = bcDesc.getSouthEastEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southFrontEdge,bcBox)){
      ierr = bcDesc.getSouthFrontEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southBackEdge,bcBox)){
      ierr = bcDesc.getSouthBackEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.westFrontEdge,bcBox)){
      ierr = bcDesc.getWestFrontEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.westBackEdge,bcBox)){
      ierr = bcDesc.getWestBackEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.eastFrontEdge,bcBox)){
      ierr = bcDesc.getEastFrontEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.eastBackEdge,bcBox)){
      ierr = bcDesc.getEastBackEdgeBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    /* Check corners */
    if (boxIntersection(tempBox,boxB.northFaceNorthWestCorner,bcBox)){
      ierr = bcDesc.getNorthFaceNorthWestCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northFaceNorthEastCorner,bcBox)){
      ierr = bcDesc.getNorthFaceNorthEastCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northFaceSouthWestCorner,bcBox)){
      ierr = bcDesc.getNorthFaceSouthWestCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.northFaceSouthEastCorner,bcBox)){
      ierr = bcDesc.getNorthFaceSouthEastCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southFaceNorthWestCorner,bcBox)){
      ierr = bcDesc.getSouthFaceNorthWestCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southFaceNorthEastCorner,bcBox)){
      ierr = bcDesc.getSouthFaceNorthEastCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southFaceSouthWestCorner,bcBox)){
      ierr = bcDesc.getSouthFaceSouthWestCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
    if (boxIntersection(tempBox,boxB.southFaceSouthEastCorner,bcBox)){
      ierr = bcDesc.getSouthFaceSouthEastCornerBoundary(bcBox,&loop); CHKERRQ(ierr);
      boundaryContainer.push_back(loop);
    }
  }
  PetscFunctionReturn(0);
}
