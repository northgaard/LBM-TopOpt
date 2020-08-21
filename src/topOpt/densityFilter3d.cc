#include "densityFilter3d.hh"

DensityFilter3d::DensityFilter3d(PetscInt filterRadius,
                                 Box3d filterDomain, const LBSolver3d& solver,
                                 const AdjointLBSolver& adjSolver)
{
  PetscErrorCode ierr;
  Box3d localDomain;
  boxIntersection(filterDomain,solver.getLocalBoundingBox(),localDomain);

  ierr = setUp(filterRadius,localDomain,
               adjSolver.sensitivityGrid);
  CHKERRABORT(solver.communicator,ierr);
}

DensityFilter3d::~DensityFilter3d()
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  if (filterMat){
    PetscObjectGetComm((PetscObject) filterMat, &comm);
    ierr = MatDestroy(&filterMat);
    CHKERRABORT(comm,ierr);
  }
  if (denominatorVec){
    PetscObjectGetComm((PetscObject) denominatorVec, &comm);
    ierr = VecDestroy(&denominatorVec);
    CHKERRABORT(comm,ierr);
  }
}

PetscErrorCode DensityFilter3d::filterSensitivities(Vec sensitivityField,
                                                    Vec* constraintSensitivityFields,
                                                    PetscInt numConstraints)
{
  PetscErrorCode ierr;
  Vec temp;
  PetscFunctionBeginUser;
  ierr = VecDuplicate(denominatorVec,&temp); CHKERRQ(ierr);
  ierr = VecPointwiseDivide(temp,sensitivityField,denominatorVec); CHKERRQ(ierr);
  ierr = MatMultTranspose(filterMat,temp,sensitivityField); CHKERRQ(ierr);
  for (PetscInt ii = 0; ii < numConstraints; ++ii){
    ierr = VecPointwiseDivide(temp,constraintSensitivityFields[ii],denominatorVec);
    CHKERRQ(ierr);
    ierr = MatMultTranspose(filterMat,temp,constraintSensitivityFields[ii]); CHKERRQ(ierr);
  }
  ierr = VecDestroy(&temp); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DensityFilter3d::setUp(PetscInt filterRadius, const Box3d& localDomain,
                                      DM theGrid)
{
  PetscErrorCode ierr;
  PetscInt Nx,Ny,Nz,NxLocal,NyLocal,NzLocal;
  DM filterGrid;
  DMBoundaryType bx,by,bz;
  DMDAStencilType sType;
  MPI_Comm comm;
  PetscFunctionBeginUser;

  PetscObjectGetComm((PetscObject) theGrid,&comm);
  ierr = DMDAGetInfo(theGrid,0,&Nx,&Ny,&Nz,&NxLocal,&NyLocal,&NzLocal,0,0,
                     &bx,&by,&bz,&sType);
  CHKERRQ(ierr);
  const PetscInt *Lx, *Ly, *Lz;
  ierr = DMDAGetOwnershipRanges(theGrid,&Lx,&Ly,&Lz); CHKERRQ(ierr);
  ierr = DMDACreate3d(comm,bx,by,bz,sType,Nx,Ny,Nz,NxLocal,NyLocal,NzLocal,
                      1,filterRadius,Lx,Ly,Lz,&filterGrid);
  CHKERRQ(ierr);

  PetscInt xs,xm,ys,ym,zs,zm;
  ierr = DMDAGetGhostCorners(filterGrid,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

  /* Create PETSc objects */
  ierr = DMSetMatType(filterGrid,MATMPIAIJ); CHKERRQ(ierr);
  ierr = DMCreateMatrix(filterGrid,&filterMat); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(filterGrid,&denominatorVec); CHKERRQ(ierr);

  /* Assemble filtering matrix */
  // Initialize as identity matrix
  Vec dummy;
  ierr = VecDuplicate(denominatorVec,&dummy); CHKERRQ(ierr);
  ierr = VecSet(dummy,1.); CHKERRQ(ierr);
  ierr = MatDiagonalSet(filterMat,dummy,INSERT_VALUES); CHKERRQ(ierr);

  // Fill in weights in filtering domain
  PetscScalar dist;
  PetscScalar floatFilterRadius = (PetscScalar) filterRadius;
  MatStencil row,col;

  PetscInt xx,yy,zz;
  PetscInt xmin = xs;
  PetscInt xmax = xs + xm - 1;
  PetscInt ymin = ys;
  PetscInt ymax = ys + ym - 1;
  PetscInt zmin = zs;
  PetscInt zmax = zs + zm - 1;

  for (auto kk : localDomain.zRange){
    for (auto jj : localDomain.yRange){
      for (auto ii : localDomain.xRange){

        row.i = ii; row.j = jj; row.k = kk;

        // We can go deeper!
        for (zz = -filterRadius; zz <= filterRadius; ++zz){
          if ((kk+zz) >= zmin && (kk+zz) <= zmax){
            for (yy = -filterRadius; yy <= filterRadius; ++yy){
              if ((jj+yy) >= ymin && (jj+yy) <= ymax){
                for (xx = -filterRadius; xx <= filterRadius; ++xx){
                  if ((ii+xx) >= xmin && (ii+xx) <= xmax){
                    col.i = ii+xx; col.j = jj+yy; col.k = kk+zz;
                    dist = PetscSqrtScalar((PetscScalar) (xx*xx + yy*yy + zz*zz));
                    if (dist <= floatFilterRadius){
                      dist = floatFilterRadius - dist;
                      MatSetValuesStencil(filterMat,1,&row,1,&col,&dist,INSERT_VALUES);
                    }
                  }
                }
              }
            }
          }
        }
        // End of crazy nesting
      }
    }
  }
  ierr = MatAssemblyBegin(filterMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(filterMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /* Sum rows of the filtering matrix for denominator */
  ierr = MatMult(filterMat,dummy,denominatorVec); CHKERRQ(ierr);

  ierr = VecDestroy(&dummy); CHKERRQ(ierr);
  ierr = DMDestroy(&filterGrid); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
