#include "fullToDesignMapping3d.hh"

FullToDesignMapping3d::FullToDesignMapping3d(Box3d designBox, LBSolver3d& solver)
  : designGrid(nullptr), worldDesignGrid(nullptr), fullToDesign(nullptr)
{
  PetscErrorCode ierr;
  /* Create design DMDA */
  communicator = solver.communicator;
  PetscInt nxd = designBox.getNx();
  PetscInt nyd = designBox.getNy();
  PetscInt nzd = designBox.getNz();
  DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE, bz = DM_BOUNDARY_NONE;
  DMDAStencilType stencilType = DMDA_STENCIL_BOX;
  PetscInt stencilWidth = 1;
  ierr = DMDACreate3d(communicator,bx,by,bz,stencilType,nxd,nyd,nzd,
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,stencilWidth,
                      0,0,0,&designGrid);
  CHKERRABORT(communicator,ierr);
  /* Create Petsc VecScatter */
  DM fullDesignGrid;
  Box3d localDesignBox;
  PetscInt lnx,lny,lnz,ltot;
  PetscInt *fullDomainIDs, *designDomainIDs;
  ierr = solver.getDesignGrid(&fullDesignGrid);
  CHKERRABORT(communicator,ierr);

  if (boxIntersection(solver.getLocalBoundingBox(),designBox,localDesignBox)){
    lnx = localDesignBox.getNx();
    lny = localDesignBox.getNy();
    lnz = localDesignBox.getNz();
  } else {
    lnx = 0;
    lny = 0;
    lnz = 0;
  }
  ltot = lnx*lny*lnz;
  ierr = PetscMalloc1(ltot,&fullDomainIDs);
  CHKERRABORT(communicator,ierr);
  ierr = PetscMalloc1(ltot,&designDomainIDs);
  CHKERRABORT(communicator,ierr);

  PetscInt shiftx = designBox.xRange.getBeginId();
  PetscInt shifty = designBox.yRange.getBeginId();
  PetscInt shiftz = designBox.zRange.getBeginId();
  PetscInt nx = solver.getBoundingBox().getNx();
  PetscInt ny = solver.getBoundingBox().getNy();
  PetscInt nScat = 0;

  for (auto kk : localDesignBox.zRange){
    for (auto jj : localDesignBox.yRange){
      for (auto ii : localDesignBox.xRange){
        fullDomainIDs[nScat] = ii + jj*nx + kk*nx*ny;
        designDomainIDs[nScat] = (ii - shiftx) + (jj - shifty)*nxd +
          (kk - shiftz)*nxd*nyd;
        ++nScat;
      }
    }
  }

  AO fullAO, designAO;
  IS fullDomainIS, designDomainIS;
  Vec fullDomainVec, designDomainVec;
  ierr = DMDAGetAO(fullDesignGrid,&fullAO); CHKERRABORT(communicator,ierr);
  ierr = DMDAGetAO(designGrid,&designAO); CHKERRABORT(communicator,ierr);
  ierr = AOApplicationToPetsc(fullAO,ltot,fullDomainIDs);
  CHKERRABORT(communicator,ierr);
  ierr = AOApplicationToPetsc(designAO,ltot,designDomainIDs);
  CHKERRABORT(communicator,ierr);
  ierr = ISCreateGeneral(communicator,ltot,fullDomainIDs,PETSC_OWN_POINTER,&fullDomainIS);
  CHKERRABORT(communicator,ierr);
  ierr = ISCreateGeneral(communicator,ltot,designDomainIDs,PETSC_OWN_POINTER,&designDomainIS);
  CHKERRABORT(communicator,ierr);
  ierr = DMCreateGlobalVector(fullDesignGrid,&fullDomainVec);
  CHKERRABORT(communicator,ierr);
  ierr = DMCreateGlobalVector(designGrid,&designDomainVec);
  CHKERRABORT(communicator,ierr);

  ierr = VecScatterCreate(fullDomainVec,fullDomainIS,designDomainVec,designDomainIS,&fullToDesign);
  CHKERRABORT(communicator,ierr);

  /* Clean up */
  ierr = ISDestroy(&fullDomainIS); CHKERRABORT(communicator,ierr);
  ierr = ISDestroy(&designDomainIS); CHKERRABORT(communicator,ierr);
  ierr = VecDestroy(&fullDomainVec); CHKERRABORT(communicator,ierr);
  ierr = VecDestroy(&designDomainVec); CHKERRABORT(communicator,ierr);
}

FullToDesignMapping3d::~FullToDesignMapping3d()
{
  PetscErrorCode ierr;
  if (designGrid){
    ierr = DMDestroy(&designGrid);
    CHKERRABORT(communicator,ierr);
  }
  if (worldDesignGrid){
    ierr = DMDestroy(&worldDesignGrid);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  if (fullToDesign){
    ierr = VecScatterDestroy(&fullToDesign);
    CHKERRABORT(communicator,ierr);
  }
}

PetscErrorCode FullToDesignMapping3d::mapFullToDesign(Vec fullDomain, Vec designDomain)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecScatterBegin(fullToDesign,fullDomain,designDomain,INSERT_VALUES,
                         SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,fullDomain,designDomain,INSERT_VALUES,
                       SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FullToDesignMapping3d::mapDesignToFull(Vec designDomain, Vec fullDomain)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecScatterBegin(fullToDesign,designDomain,fullDomain,INSERT_VALUES,
                         SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,designDomain,fullDomain,INSERT_VALUES,
                       SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FullToDesignMapping3d::createDesignDomainVec(Vec* vec)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = DMCreateGlobalVector(designGrid,vec); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  This creates a design domain vector on PETSC_COMM_WORLD, regardless of the communicator
  otherwise used. If the communicator is indeed PETSC_COMM_WORLD, this just calls
  createDesignDomainVec()
*/
PetscErrorCode FullToDesignMapping3d::createWorldDesignDomainVec(Vec* vec)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  if (communicator == PETSC_COMM_WORLD){
    ierr = createDesignDomainVec(vec); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  // World design DMDA is lazily initialized
  if (!worldDesignGrid){
    PetscInt nx,ny,nz,dof,stencilWidth;
    DMBoundaryType bx,by,bz;
    DMDAStencilType stencilType;
    ierr = DMDAGetInfo(designGrid,0,&nx,&ny,&nz,0,0,0,&dof,&stencilWidth,&bx,&by,
                       &bz,&stencilType); CHKERRQ(ierr);
    ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stencilType,nx,ny,nz,
                        PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof,stencilWidth,
                        0,0,0,&worldDesignGrid); CHKERRQ(ierr);
  }
  ierr = DMCreateGlobalVector(worldDesignGrid,vec); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
