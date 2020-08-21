#include "densityFilter2d.hh"

DensityFilter2d::DensityFilter2d(PetscInt filterRadius,
                                 Box2d filterDomain, const NewLBSolver2d& solver,
                                 const AdjointLBSolver& adjSolver)

{
  PetscErrorCode ierr;
  Box2d localDomain;
  boxIntersection(filterDomain,solver.getLocalBoundingBox(),localDomain);

  ierr = setUp(filterRadius,localDomain,
	       adjSolver.sensitivityGrid);
  CHKERRABORT(solver.communicator,ierr);
}

DensityFilter2d::~DensityFilter2d()
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

PetscErrorCode DensityFilter2d::filterDesign(Vec designField, Vec filteredField)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = MatMult(filterMat,designField,filteredField); CHKERRQ(ierr);
  ierr = VecPointwiseDivide(filteredField,filteredField,denominatorVec);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DensityFilter2d::filterSensitivities(Vec sensitivityField,
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

PetscErrorCode DensityFilter2d::setUp(PetscInt filterRadius,
                                      const Box2d& localDomain,
                                      DM theGrid)
{
  PetscErrorCode ierr;
  PetscInt Nx,Ny,NxLocal,NyLocal;
  DM filterGrid;
  DMBoundaryType bx,by;
  DMDAStencilType sType;
  MPI_Comm comm;
  PetscFunctionBeginUser;

  PetscObjectGetComm((PetscObject) theGrid, &comm);
  ierr = DMDAGetInfo(theGrid,0,&Nx,&Ny,0,&NxLocal,&NyLocal,0,0,0,&bx,&by,
		     0,&sType);
  CHKERRQ(ierr);
  const PetscInt* Lx, *Ly;
  ierr = DMDAGetOwnershipRanges(theGrid,&Lx,&Ly,0); CHKERRQ(ierr);
  ierr = DMDACreate2d(comm,bx,by,sType,Nx,Ny,NxLocal,NyLocal,
		      1,filterRadius,Lx,Ly,&filterGrid);
  CHKERRQ(ierr);

  PetscInt xs,xm,ys,ym;
  ierr = DMDAGetGhostCorners(filterGrid,&xs,&ys,0,&xm,&ym,0); CHKERRQ(ierr);

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

  PetscInt xx,yy;
  PetscInt xmin = xs;
  PetscInt xmax = xs + xm - 1;
  PetscInt ymin = ys;
  PetscInt ymax = ys + ym - 1;

  for (auto jj : localDomain.yRange){
    for (auto ii : localDomain.xRange){

      row.i = ii; row.j = jj;

      for (yy = -filterRadius; yy <= filterRadius; ++yy){
        if ((jj+yy) >= ymin && (jj+yy) <= ymax){
          for (xx = -filterRadius; xx <= filterRadius; ++xx){
            if ((ii+xx) >= xmin && (ii+xx) <= xmax){
              col.i = ii+xx; col.j = jj+yy;

              dist = PetscSqrtScalar((PetscScalar) (xx*xx + yy*yy));

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
  ierr = MatAssemblyBegin(filterMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(filterMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /* Sum rows of the filtering matrix for denominator */
  ierr = MatMult(filterMat,dummy,denominatorVec); CHKERRQ(ierr);

  ierr = VecDestroy(&dummy); CHKERRQ(ierr);
  ierr = DMDestroy(&filterGrid); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
