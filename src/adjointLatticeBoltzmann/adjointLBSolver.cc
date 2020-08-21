#include "adjointLBSolver.hh"

AdjointLBSolver::AdjointLBSolver() : communicator(MPI_COMM_NULL),
                                     adjointGlobal(nullptr),
                                     sensitivityGlobal(nullptr), latticeGrid(nullptr),
                                     materialGrid(nullptr), sensitivityGrid(nullptr),
                                     adjointCollideLoop(nullptr), adjointStreamLoop(nullptr)
{}

AdjointLBSolver::~AdjointLBSolver()
{
  PetscErrorCode ierr;
  if (adjointGlobal){
    ierr = VecDestroy(&adjointGlobal);
    CHKERRABORT(communicator,ierr);
  }
  if (sensitivityGlobal){
    ierr = VecDestroy(&sensitivityGlobal);
    CHKERRABORT(communicator,ierr);
  }
  if (adjointCollideLoop){ delete adjointCollideLoop; }
  if (adjointStreamLoop){ delete adjointStreamLoop; }
  for (const auto& adjointBoundary : adjointBoundaryContainer){
    if (adjointBoundary){ delete adjointBoundary; }
  }
}

PetscErrorCode AdjointLBSolver::adjointCollide(Vec distributions,
                                               Vec obstacle)
{
  PetscErrorCode ierr;
  void *adj, *fdist, *obst, *sens;
  PetscFunctionBeginUser;

  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(materialGrid,obstacle,&obst); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(sensitivityGrid,sensitivityGlobal,&sens); CHKERRQ(ierr);

  adjointCollideLoop->execute(adj,fdist,obst,sens);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(materialGrid,obstacle,&obst); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(sensitivityGrid,sensitivityGlobal,&sens); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode AdjointLBSolver::adjointStream(Vec distributions)
{
  PetscErrorCode ierr;
  void *adj, *adjPrev, *fdist;
  PetscFunctionBeginUser;

  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);

  adjointBoundaryComputations(adj,fdist);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);

  Vec adjointLocal;
  ierr = DMGetLocalVector(latticeGrid,&adjointLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(latticeGrid,adjointGlobal,INSERT_VALUES,
                              adjointLocal); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(latticeGrid,adjointGlobal,INSERT_VALUES,
                            adjointLocal); CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointLocal,&adjPrev); CHKERRQ(ierr);

  adjointStreamLoop->execute(adj,adjPrev,nullptr,nullptr);
  --currentTimestep;

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointLocal,&adjPrev); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(latticeGrid,&adjointLocal); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode AdjointLBSolver::computeCollideSource(const LBMacroFunctional& func,
                                                     const PetscScalar timeConst,
                                                     Vec distributions,
                                                     Vec obstacle)
{
  PetscErrorCode ierr;
  void *adj, *fdist, *obst;
  PetscFunctionBeginUser;
  ierr = DMDAVecGetArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(materialGrid,obstacle,&obst); CHKERRQ(ierr);

  func.adjointCollideSource(currentTimestep,timeConst,adj,fdist,obst);

  ierr = DMDAVecRestoreArrayDOF(latticeGrid,adjointGlobal,&adj); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(materialGrid,obstacle,&obst); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode AdjointLBSolver::computeSensitivitySource(const LBMacroFunctional& func,
                                                         const PetscScalar timeConst,
                                                         Vec distributions,
                                                         Vec obstacle)
{
  PetscErrorCode ierr;
  void *sens, *fdist, *obst;
  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(sensitivityGrid,sensitivityGlobal,&sens); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(materialGrid,obstacle,&obst); CHKERRQ(ierr);

  func.adjointSensitivitySource(currentTimestep,timeConst,sens,fdist,obst);

  ierr = DMDAVecRestoreArray(sensitivityGrid,sensitivityGlobal,&sens); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(latticeGrid,distributions,&fdist); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(materialGrid,obstacle,&obst); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode AdjointLBSolver::outputAdjoints(const std::string& folder,
                                               const std::string& name)
{
  PetscErrorCode ierr;
  PetscViewer view;
  PetscFunctionBeginUser;
  std::string output;
  output = folder + name;
  ierr = PetscViewerVTKOpen(communicator,output.c_str(),FILE_MODE_WRITE,&view); CHKERRQ(ierr);
  ierr = VecView(adjointGlobal,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode AdjointLBSolver::resetAdjoints()
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecSet(adjointGlobal,0.); CHKERRQ(ierr);
  ierr = VecSet(sensitivityGlobal,0.); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode AdjointLBSolver::setCurrentTimestep(PetscInt _ct)
{
  currentTimestep = _ct;
  return 0;
}

PetscErrorCode AdjointLBSolver::getSensitivities(Vec* theVec)
{
  *theVec = sensitivityGlobal;
  return 0;
}

PetscErrorCode AdjointLBSolver::createPetscObjects(DM _lat, DM _mat, DM _sens)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  latticeGrid = _lat;
  materialGrid = _mat;
  sensitivityGrid = _sens;
  ierr = DMCreateGlobalVector(latticeGrid,&adjointGlobal); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(sensitivityGrid,&sensitivityGlobal); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void AdjointLBSolver::adjointBoundaryComputations(void* adj, void* fdist)
{
  for (const auto& boundary : adjointBoundaryContainer){
    boundary->execute(currentTimestep,adj,fdist);
  }
}
