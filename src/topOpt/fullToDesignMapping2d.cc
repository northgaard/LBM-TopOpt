#include "fullToDesignMapping2d.hh"

FullToDesignMapping2d::FullToDesignMapping2d(Box2d designBox, NewLBSolver2d& solver)
  : designGrid(nullptr), worldDesignGrid(nullptr), fullToDesign(nullptr)
{
  PetscErrorCode ierr;
  /* Create Petsc design DMDA */
  communicator = solver.communicator;
  PetscInt nxd = designBox.getNx();
  PetscInt nyd = designBox.getNy();
  DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE;
  DMDAStencilType stencilType = DMDA_STENCIL_BOX;
  PetscInt stencilWidth = 1;
  ierr = DMDACreate2d(communicator,bx,by,stencilType,nxd,nyd,PETSC_DECIDE,
                      PETSC_DECIDE,1,stencilWidth,0,0,&designGrid);
  CHKERRABORT(communicator,ierr);
  /* Create Petsc VecScatter */
  DM fullDesignGrid;
  Box2d localDesignBox;
  PetscInt lnx,lny,ltot;
  PetscInt *fullDomainIDs, *designDomainIDs;
  ierr = solver.getDesignGrid(&fullDesignGrid);
  CHKERRABORT(communicator,ierr);

  if (boxIntersection(solver.getLocalBoundingBox(),designBox,localDesignBox)){
    lnx = localDesignBox.getNx();
    lny = localDesignBox.getNy();
  } else {
    lnx = 0;
    lny = 0;
  }
  ltot = lnx*lny;
  ierr = PetscMalloc1(ltot,&fullDomainIDs);
  CHKERRABORT(communicator,ierr);
  ierr = PetscMalloc1(ltot,&designDomainIDs);
  CHKERRABORT(communicator,ierr);

  PetscInt shiftx = designBox.xRange.getBeginId();
  PetscInt shifty = designBox.yRange.getBeginId();
  PetscInt nx = solver.getBoundingBox().getNx();
  PetscInt nScat = 0;

  for (auto jj : localDesignBox.yRange){
    for (auto ii : localDesignBox.xRange){
      fullDomainIDs[nScat] = ii + jj*nx;
      designDomainIDs[nScat] = (ii - shiftx) + (jj - shifty)*nxd;
      ++nScat;
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

FullToDesignMapping2d::~FullToDesignMapping2d()
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

PetscErrorCode FullToDesignMapping2d::mapFullToDesign(Vec fullDomain, Vec designDomain)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecScatterBegin(fullToDesign,fullDomain,designDomain,INSERT_VALUES,
                         SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,fullDomain,designDomain,INSERT_VALUES,
                       SCATTER_FORWARD); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FullToDesignMapping2d::mapDesignToFull(Vec designDomain, Vec fullDomain)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecScatterBegin(fullToDesign,designDomain,fullDomain,INSERT_VALUES,
                         SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(fullToDesign,designDomain,fullDomain,INSERT_VALUES,
                       SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FullToDesignMapping2d::createDesignDomainVec(Vec* vec)
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
PetscErrorCode FullToDesignMapping2d::createWorldDesignDomainVec(Vec* vec)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  if (communicator == PETSC_COMM_WORLD){
    ierr = createDesignDomainVec(vec); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  // World design DMDA is lazily initialized
  if (!worldDesignGrid){
    PetscInt nx,ny,dof,stencilWidth;
    DMBoundaryType bx,by;
    DMDAStencilType stencilType;
    ierr = DMDAGetInfo(designGrid,0,&nx,&ny,0,0,0,0,&dof,&stencilWidth,&bx,&by,
                       0,&stencilType); CHKERRQ(ierr);
    ierr = DMDACreate2d(PETSC_COMM_WORLD,bx,by,stencilType,nx,ny,PETSC_DECIDE,
                        PETSC_DECIDE,dof,stencilWidth,0,0,&worldDesignGrid); CHKERRQ(ierr);
  }
  ierr = DMCreateGlobalVector(worldDesignGrid,vec); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
