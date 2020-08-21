#include "LBSolver.hh"

LBSolver::LBSolver() : communicator(MPI_COMM_NULL),
                       distributionsGlobal(nullptr), macrosGlobal(nullptr),
                       initMacrosGlobal(nullptr), materialGlobal(nullptr),
                       latticeGrid(nullptr), macroGrid(nullptr), materialGrid(nullptr),
                       collideAndSwapLoop(nullptr),collideAndStreamLoop(nullptr),
                       collideLoop(nullptr),streamBySwappingLoop(nullptr),
                       streamLoop(nullptr),computeMacrosLoop(nullptr),
                       initializeAtEquilibriumLoop(nullptr),
                       currentState(INVALID), currentTimestep(0),
                       designGrid(nullptr)
{}

LBSolver::~LBSolver()
{
  PetscErrorCode ierr;
  if (distributionsGlobal){
    ierr = VecDestroy(&distributionsGlobal);
    CHKERRABORT(communicator,ierr);
  }
  if (macrosGlobal){
    ierr = VecDestroy(&macrosGlobal);
    CHKERRABORT(communicator,ierr);
  }
  if (initMacrosGlobal){
    ierr = VecDestroy(&initMacrosGlobal);
    CHKERRABORT(communicator,ierr);
  }
  if (materialGlobal){
    ierr = VecDestroy(&materialGlobal);
    CHKERRABORT(communicator,ierr);
  }
  if (latticeGrid){
    ierr = DMDestroy(&latticeGrid);
    CHKERRABORT(communicator,ierr);
  }
  if (macroGrid){
    ierr = DMDestroy(&macroGrid);
    CHKERRABORT(communicator,ierr);
  }
  if (materialGrid){
    if (materialGrid == designGrid){
      designGrid = nullptr;
    }
    ierr = DMDestroy(&materialGrid);
    CHKERRABORT(communicator,ierr);
  }
  if (designGrid){
    ierr = DMDestroy(&designGrid);
    CHKERRABORT(communicator,ierr);
  }
  if (collideAndSwapLoop){ delete collideAndSwapLoop; }
  if (collideAndStreamLoop){ delete collideAndStreamLoop; }
  if (collideLoop){ delete collideLoop; }
  if (streamBySwappingLoop){ delete streamBySwappingLoop; }
  if (streamLoop){ delete streamLoop; }
  if (computeMacrosLoop){ delete computeMacrosLoop; }
  if (initializeAtEquilibriumLoop){ delete initializeAtEquilibriumLoop; }
  for (const auto& boundary : boundaryContainer){
    if (boundary){ delete boundary; }
  }
}

PetscErrorCode LBSolver::setDistributions(Vec newDist)
{
  if (newDist == distributionsGlobal){ return 0; }
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecCopy(newDist,distributionsGlobal); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::computeMacroFunctional(const LBMacroFunctional& func,
                                                const PetscScalar timeConst,
                                                PetscScalar* value)
{
  PetscErrorCode ierr;
  void *mac, *obst;
  PetscFunctionBeginUser;
  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  if (materialGrid){
    ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&obst); CHKERRQ(ierr);
  }
  func.evaluate(currentTimestep,timeConst,mac,obst,value);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&mac); CHKERRQ(ierr);
  if (materialGrid){
    ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&obst); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::setCurrentTimestep(PetscInt _ct)
{
  currentTimestep = _ct;
  return 0;
}

PetscErrorCode LBSolver::streamAndCollide()
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = stream(); CHKERRQ(ierr);
  ierr = collide(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::collideAndStream()
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = collide(); CHKERRQ(ierr);
  ierr = stream(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// PetscErrorCode LBSolver::collideAndStream()
// {
//   PetscErrorCode ierr;
//   void *distPtr, *macroPtr, *obstaclePtr;
//   ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
//   ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
//   if (materialGrid){
//     ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
//   }

//   collideAndStreamLoop->execute(distPtr,macroPtr,obstaclePtr);
//   ++currentTimestep;
//   boundaryComputations(distPtr);

//   ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
//   ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
//   if (materialGrid){
//     ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
//   }
//   ierr = DMLocalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
//                              distributionsGlobal); CHKERRQ(ierr);
//   ierr = DMLocalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
//                            distributionsGlobal); CHKERRQ(ierr);
//   return 0;
// }

// PetscErrorCode LBSolver::collideAndSwap()
// {
//   PetscErrorCode ierr;
//   void *distPtr, *macroPtr, *obstaclePtr;
//   ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
//   ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
//   if (materialGrid){
//     ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
//   }

//   collideAndSwapLoop->execute(distPtr,macroPtr,obstaclePtr);

//   ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
//   ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
//   if (materialGrid){
//     ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
//   }
//   return 0;
// }

PetscErrorCode LBSolver::collide()
{
  PetscErrorCode ierr;
  void *distPtr, *macroPtr, *obstaclePtr;
  PetscFunctionBeginUser;
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
  if (materialGrid){
    ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
  }

  collideLoop->execute(distPtr,macroPtr,obstaclePtr);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
  if (materialGrid){
    ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::stream()
{
  PetscErrorCode ierr;
  void *distPtr, *localDistPtr;
  Vec distributionsLocal;
  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(latticeGrid,&distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
                              distributionsLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
                            distributionsLocal); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsLocal,&localDistPtr); CHKERRQ(ierr);

  streamLoop->execute(distPtr,localDistPtr,nullptr);
  ++currentTimestep;
  boundaryComputations(distPtr);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsLocal,&localDistPtr); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(latticeGrid,&distributionsLocal); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// PetscErrorCode LBSolver::streamBySwapping()
// {
//   PetscErrorCode ierr;
//   void *distPtr, *macroPtr, *obstaclePtr;
//   ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
//   ierr = DMDAVecGetArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
//   if (materialGrid){
//     ierr = DMDAVecGetArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
//   }

//   streamBySwappingLoop->execute(distPtr,nullptr,nullptr);
//   ++currentTimestep;
//   boundaryComputations(distPtr);

//   ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
//   ierr = DMDAVecRestoreArrayDOF(macroGrid,macrosGlobal,&macroPtr); CHKERRQ(ierr);
//   if (materialGrid){
//     ierr = DMDAVecRestoreArrayDOF(materialGrid,materialGlobal,&obstaclePtr); CHKERRQ(ierr);
//   }
//   ierr = DMLocalToLocalBegin(latticeGrid,distributionsGlobal,INSERT_VALUES,
//                              distributionsGlobal); CHKERRQ(ierr);
//   ierr = DMLocalToLocalEnd(latticeGrid,distributionsGlobal,INSERT_VALUES,
//                            distributionsGlobal); CHKERRQ(ierr);
//   return 0;
// }

PetscErrorCode LBSolver::initializeAtEquilibrium()
{
  PetscErrorCode ierr;
  void *distPtr, *initMacroPtr;

  ierr = DMDAVecGetArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(macroGrid,initMacrosGlobal,&initMacroPtr); CHKERRQ(ierr);

  initializeAtEquilibriumLoop->execute(distPtr,initMacroPtr,nullptr);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributionsGlobal,&distPtr); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(macroGrid,initMacrosGlobal,&initMacroPtr); CHKERRQ(ierr);
  ierr = VecCopy(initMacrosGlobal,macrosGlobal); CHKERRQ(ierr);
  currentTimestep = 0;
  return 0;
}

PetscErrorCode LBSolver::outputMacros(const std::string& folder, const std::string& name)
{
  PetscErrorCode ierr;
  PetscViewer view;
  const std::string output = folder + name;
  PetscFunctionBeginUser;
  ierr = PetscViewerVTKOpen(communicator,output.c_str(),FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(macrosGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::outputDistributions(const std::string& folder, const std::string& name)
{
  PetscErrorCode ierr;
  PetscViewer view;
  const std::string output = folder + name;
  PetscFunctionBeginUser;
  ierr = PetscViewerVTKOpen(communicator,output.c_str(),FILE_MODE_WRITE,&view);
  CHKERRQ(ierr);
  ierr = VecView(distributionsGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::outputFields(const std::string& folder, const std::string& name)
{
  PetscErrorCode ierr;
  PetscViewer view;
  const std::string output = folder + name;
  ierr = PetscViewerVTKOpen(communicator,output.c_str(),FILE_MODE_WRITE,&view); CHKERRQ(ierr);
  ierr = VecView(materialGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void LBSolver::boundaryComputations(void* distPtr)
{
  for (const auto& boundary : boundaryContainer){
    boundary->execute(currentTimestep,distPtr,distPtr);
  }
}

PetscErrorCode LBSolver::getDesignGrid(DM* dm)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  if (!designGrid){ ierr = initializeDesignGrid(); CHKERRQ(ierr); }
  if (dm){ *dm = designGrid; }
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::initializeDesignGrid()
{
  PetscErrorCode ierr;
  PetscInt materialDOF;
  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(materialGrid,0,0,0,0,0,0,0,&materialDOF,0,0,0,0,0); CHKERRQ(ierr);
  if (materialDOF == 1){
    designGrid = materialGrid;
  } else {
    ierr = DMDAGetReducedDMDA(latticeGrid,1,&designGrid); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::createDesignVec(Vec* theVec)
{
  PetscErrorCode ierr;
  DM dg;
  PetscFunctionBeginUser;
  ierr = getDesignGrid(&dg); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dg,theVec); CHKERRQ(ierr);
  if (materialGrid == dg){
    ierr = VecCopy(materialGlobal,*theVec); CHKERRQ(ierr);
  } else {
    // TODO: implement this case
  }
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::setFieldFromDesignVec(Vec theVec)
{
  PetscErrorCode ierr;
  DM dg;
  PetscFunctionBeginUser;
  ierr = getDesignGrid(&dg); CHKERRQ(ierr);
  // TODO: error check on layout of theVec
  if (materialGrid == dg){
    ierr = VecCopy(theVec,materialGlobal); CHKERRQ(ierr);
  } else {
    PetscInt idDesign = 0;
    ierr = VecStrideSubSetScatter(theVec,1,&idDesign,&idDesign,materialGlobal,INSERT_VALUES);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode LBSolver::setUpAdjointSolver(AdjointLBSolver* adjSolver)
{
  PetscErrorCode ierr;
  adjSolver->communicator = communicator;
  ierr = adjSolver->createPetscObjects(latticeGrid,materialGrid,designGrid); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode LBSolver::setAdjointCollisionLoop(AdjointLBLoop* adjLoop,
                                                 AdjointLBSolver* adjSolver)
{
  adjSolver->adjointCollideLoop = adjLoop;
  return 0;
}

PetscErrorCode LBSolver::setAdjointStreamingLoop(AdjointLBLoop* adjLoop,
                                              AdjointLBSolver* adjSolver)
{
  adjSolver->adjointStreamLoop = adjLoop;
  return 0;
}

PetscErrorCode LBSolver::addAdjointBoundaryCondition(AdjointLBBoundaryLoop* adjLoop,
                                                     AdjointLBSolver* adjSolver)
{
  adjSolver->adjointBoundaryContainer.push_back(adjLoop);
  return 0;
}
