#include "LBTopOpt2d.hh"

LBTopOpt2d::LBTopOpt2d(ObstacleLBSolver2d& _solv, AdjointLBSolver2d& _asolv,
		       std::unique_ptr<Filter> _fil)
  : solver(_solv), adjointSolver(_asolv), filter(std::move(_fil))
{

  PetscErrorCode ierr;
  PetscInt obstacleDOF;
  ierr = DMDAGetInfo(solver.obstacleGrid,0,0,0,0,0,0,0,&obstacleDOF,0,0,0,0,0);

  if (obstacleDOF == 1){
    ierr = VecDuplicate(solver.obstacleGlobal,&fullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) fullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDuplicate(solver.obstacleGlobal,&physicalFullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) physicalFullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDuplicate(solver.obstacleGlobal,&sensitivityFull);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) sensitivityFull);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);

    ierr = VecCopy(solver.obstacleGlobal,fullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecCopy(solver.obstacleGlobal,physicalFullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);

    
  } else {
    DM singleGrid;
    ierr = DMDAGetReducedDMDA(solver.obstacleGrid,1,&singleGrid);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = DMCreateGlobalVector(singleGrid,&fullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) fullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = DMCreateGlobalVector(singleGrid,&physicalFullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) physicalFullDomainVec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = DMCreateGlobalVector(singleGrid,&sensitivityFull);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) sensitivityFull);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = DMDestroy(&singleGrid);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }

  constraintValues.push_back(0.);

}

PetscErrorCode LBTopOpt2d::
getObjectiveAndConstraints(PetscScalar& obj, PetscScalar& con)
{

  PetscErrorCode ierr;
  obj = objectiveValue;

  ierr = VecScatterBegin(fullToDesign,physicalFullDomainVec,physicalDesignDomainVec,
			 INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,physicalFullDomainVec,physicalDesignDomainVec,
		       INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = constraint.computeConstraint(physicalDesignDomainVec,con); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBTopOpt2d::
getObjectiveConstraintsAndSensitivities(PetscScalar& obj, PetscScalar& con,
					Vec& objSens, std::vector<Vec>& consSens)
{

  PetscErrorCode ierr;

  ierr = getObjectiveAndConstraints(obj,con); CHKERRQ(ierr);

  ierr = VecScatterBegin(fullToDesign,sensitivityFull,sensitivityDesign,
			 INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,sensitivityFull,sensitivityDesign,
		       INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  objSens = sensitivityDesign;

  ierr = constraint.computeSensitivities(constraintSensitivities[0]);
  CHKERRQ(ierr);

  if (filter){

    Vec temp;
    ierr = VecDuplicate(sensitivityFull,&temp); CHKERRQ(ierr);
    ierr = VecSet(temp,0.); CHKERRQ(ierr);

    ierr = VecScatterBegin(fullToDesign,constraintSensitivities[0],temp,
			   INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(fullToDesign,constraintSensitivities[0],temp,
			 INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    
    ierr = filter->filterSensitivities(temp); CHKERRQ(ierr);

    ierr = VecScatterBegin(fullToDesign,temp,constraintSensitivities[0],
			   INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(fullToDesign,temp,constraintSensitivities[0],
			 INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecDestroy(&temp); CHKERRQ(ierr);
    
  }

  consSens = constraintSensitivities;
  
  return 0;

}

PetscErrorCode LBTopOpt2d::setDesignDomain(const Box2d designBox)
{

  PetscErrorCode ierr;
  ierr = createDesignObjects(designBox); CHKERRQ(ierr);
  DM singleGrid;
  AO fullAO, designAO;
  IS fullDomainIS, designDomainIS;
  PetscInt* fullDomainIDs;
  PetscInt* designDomainIDs;
  PetscInt numx,numy,numTot,obstacleDOF;
  Box2d localDesignBox;

  ierr = DMDAGetInfo(solver.obstacleGrid,0,0,0,0,0,0,0,&obstacleDOF,0,0,0,0,0);
  CHKERRQ(ierr);

  if (obstacleDOF == 1){
    singleGrid = solver.obstacleGrid;
  } else {
    ierr = DMDAGetReducedDMDA(solver.obstacleGrid,1,&singleGrid); CHKERRQ(ierr);
  }

  if (boxIntersection(solver.getLocalBoundingBox(),designBox,localDesignBox)){
    numx = localDesignBox.getNx();
    numy = localDesignBox.getNy();
  } else {
    numx = 0;
    numy = 0;
  }
  numTot = numx*numy;

  ierr = PetscMalloc1(numTot,&fullDomainIDs); CHKERRQ(ierr);
  ierr = PetscMalloc1(numTot,&designDomainIDs); CHKERRQ(ierr);

  PetscInt shiftx = designBox.xRange.getBeginId();
  PetscInt shifty = designBox.yRange.getBeginId();
  PetscInt nx = solver.getBoundingBox().getNx();
  PetscInt nxd = designBox.getNx();

  PetscInt nScat = 0;
  for (auto jj : localDesignBox.yRange){
    for (auto ii : localDesignBox.xRange){
      fullDomainIDs[nScat] = ii + jj*nx;
      designDomainIDs[nScat] = (ii - shiftx) + (jj - shifty)*nxd;
      ++nScat;
    }
  }

  ierr = DMDAGetAO(singleGrid,&fullAO); CHKERRQ(ierr);
  ierr = DMDAGetAO(designGrid,&designAO); CHKERRQ(ierr);
  ierr = AOApplicationToPetsc(fullAO,numTot,fullDomainIDs); CHKERRQ(ierr);
  ierr = AOApplicationToPetsc(designAO,numTot,designDomainIDs); CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,numTot,fullDomainIDs,
			 PETSC_COPY_VALUES,&fullDomainIS);
  CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD,numTot,designDomainIDs,
			 PETSC_COPY_VALUES,&designDomainIS);
  CHKERRQ(ierr);

  ierr = VecScatterCreate(fullDomainVec,fullDomainIS,
			  designDomainVec,designDomainIS,&fullToDesign);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) fullToDesign);
  CHKERRQ(ierr);

  ierr = PetscFree(fullDomainIDs); CHKERRQ(ierr);
  ierr = PetscFree(designDomainIDs); CHKERRQ(ierr);
  ierr = ISDestroy(&fullDomainIS); CHKERRQ(ierr);
  ierr = ISDestroy(&designDomainIS); CHKERRQ(ierr);

  if (obstacleDOF != 1){
    ierr = DMDestroy(&singleGrid); CHKERRQ(ierr);
  }

  return 0;

}

PetscErrorCode LBTopOpt2d::createDesignObjects(const Box2d& designBox)
{

  PetscErrorCode ierr;
  PetscInt nx = designBox.getNx();
  PetscInt ny = designBox.getNy();

  #ifdef PETSC35
  DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE;
  #else
  DMDABoundaryType bx = DMDA_BOUNDARY_NONE, by = DMDA_BOUNDARY_NONE;
  #endif
  DMDAStencilType stencilType = DMDA_STENCIL_BOX;
  PetscInt stencilWidth = 1;

  ierr = DMDACreate2d(PETSC_COMM_WORLD,bx,by,stencilType,nx,ny,PETSC_DECIDE,
		      PETSC_DECIDE,1,stencilWidth,0,0,&designGrid);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) designGrid);
  CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(designGrid,&designDomainVec);
  CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) designDomainVec);
  CHKERRQ(ierr);

  ierr = VecDuplicate(designDomainVec,&sensitivityDesign); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) sensitivityDesign);
  CHKERRQ(ierr);

  ierr = VecDuplicate(designDomainVec,&physicalDesignDomainVec); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) physicalDesignDomainVec);
  CHKERRQ(ierr);

  Vec temp;
  ierr = VecDuplicate(designDomainVec,&temp); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) temp); CHKERRQ(ierr);
  constraintSensitivities.push_back(temp);

  return 0;

}

PetscErrorCode LBTopOpt2d::getDesignVector(Vec& theVec)
{

  PetscErrorCode ierr;
  
  ierr = VecScatterBegin(fullToDesign,fullDomainVec,designDomainVec,INSERT_VALUES,
			 SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,fullDomainVec,designDomainVec,INSERT_VALUES,
		       SCATTER_FORWARD); CHKERRQ(ierr);

  theVec = designDomainVec;

  return 0;

}

PetscErrorCode LBTopOpt2d::restoreDesignVector(const Vec& theVec)
{

  PetscErrorCode ierr;

  ierr = VecScatterBegin(fullToDesign,theVec,fullDomainVec,INSERT_VALUES,
			 SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,theVec,fullDomainVec,INSERT_VALUES,
		       SCATTER_REVERSE); CHKERRQ(ierr);
  
  return 0;

}

PetscErrorCode LBTopOpt2d::forwardAllocation(PetscInt numCheckPoints,
					     ForwardSolutionStorage* store)
{

  PetscErrorCode ierr;
  ierr = store->allocate(numCheckPoints,solver.distributionsGlobal); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBTopOpt2d::resetSensitivities()
{

  PetscErrorCode ierr;
  ierr = VecSet(sensitivityFull,0.); CHKERRQ(ierr);
  ierr = VecSet(sensitivityDesign,0.); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode LBTopOpt2d::
addVolumeConstraint(PetscScalar _vol, VolumeConstraintType _t)
{
  constraint = VolumeConstraint(_vol,_t);

  PetscErrorCode ierr;
  Vec temp;
  ierr = DMCreateGlobalVector(designGrid,&temp); CHKERRQ(ierr);
  ierr = PetscObjectRegisterDestroy((PetscObject) temp); CHKERRQ(ierr);
  constraintSensitivities.push_back(temp);
    
  return 0;
}

std::unique_ptr<MMA> LBTopOpt2d::getMMA()
{

  PetscInt numDesign, numConstraints;
  VecGetSize(designDomainVec,&numDesign);
  numConstraints = constraintSensitivities.size();

  return std::unique_ptr<MMA>(new MMA(numDesign,numConstraints,designDomainVec));

}

PetscErrorCode LBTopOpt2d::outputFullDomain(PetscInt iter)
{

  PetscErrorCode ierr;
  PetscViewer view;
  char output[50];

  if (filter){
    ierr = filter->filterDesign(fullDomainVec,physicalFullDomainVec);
    CHKERRQ(ierr);
  }
  
  sprintf(output,"fullDomain_iter_%i.vts",iter);

  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,output,FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(physicalFullDomainVec,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);

  return 0;

}

