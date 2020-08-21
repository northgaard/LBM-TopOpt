#include "TOpTEnLB3d.hh"

static char help[] = "Implementation of the diffuser problem in 3d.";
static char outputFolder[] = "/home/sebnorg/WorkData/Diffuser3d/";

using Lattice = D3Q19;
// using CollisionOperator = StandardBGK3d<Lattice>;
using CollisionOperator = StandardMRT3d<Lattice>;
// using CollisionOperator = IncompressibleMRT3d<Lattice>;

struct VelocityOp {
  void operator()(PetscInt,PetscInt,PetscInt,PetscInt,
                  PetscScalar& ux, PetscScalar& uy, PetscScalar& uz) const
  {
    ux = vmax;
    uy = 0.;
    uz = 0.;
  }
  PetscScalar vmax;
};

struct DevFlowYZ {
  DevFlowYZ(const Box3d& area, PetscScalar _vm) : vmax(_vm)
  {
    ymid = 0.5*(area.yRange.getEndId() - area.yRange.getBeginId());
    zmid = 0.5*(area.zRange.getEndId() - area.zRange.getBeginId());
    rad = 0.5*(area.getNy() - 1);
  }
  void operator()(PetscInt,PetscInt,PetscInt idy,PetscInt idz,
                  PetscScalar& ux, PetscScalar& uy, PetscScalar& uz) const
  {
    PetscScalar circ = (idy - ymid)*(idy - ymid) + (idz - zmid)*(idz - zmid);
    if (circ <= rad*rad){
      ux = vmax*(1. - (1./(rad*rad))*circ);
    } else {
      ux = 0.;
    }
    uy = 0.;
    uz = 0.;
  }
private:
  PetscScalar vmax;
  PetscScalar ymid;
  PetscScalar zmid;
  PetscScalar rad;
};

struct PressureOp {
  void operator()(PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar& rho) const
  {
    rho = rho0;
  }
  PetscScalar rho0;
};

int main(int argc, char *argv[])
{
  TopTenInit init(argc,argv,help);
  PetscErrorCode ierr;
  IsothermalLBSolver3d<CollisionOperator> solver;

  LBSolverInfo3d info;
  PetscInt nx = 126;
  PetscInt ny = 42;
  PetscInt nz = 42;
  info.nx = nx;
  info.ny = ny;
  info.nz = nz;
  PetscInt th = ny / 3;
  IncompressibleFlowParameters par;
  par.velocityChar = 0.1;
  par.ReynoldsNumber = 100;
  par.lengthChar = ny;
  CollisionOperator op(par);
  solver.make(PETSC_COMM_WORLD,info,op);

  /* Define boundaries */
  OnGridBounceback3d<Lattice> bb;
  ierr = solver.addBoundaryCondition(Box3d(0,nx-1,0,0,0,nz-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(0,nx-1,ny-1,ny-1,0,nz-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(0,nx-1,1,ny-2,0,0),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(0,nx-1,1,ny-2,nz-1,nz-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(nx-1,nx-1,1,th,1,nz-2),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(nx-1,nx-1,2*th+1,ny-2,1,nz-2),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(nx-1,nx-1,th+1,2*th,1,th),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(nx-1,nx-1,th+1,2*th,2*th+1,nz-2),bb); CHKERRQ(ierr);

  // VelocityOp vel;
  // vel.vmax = par.velocityChar;
  DevFlowYZ vel(Box3d(0,0,1,ny-2,1,nz-2),par.velocityChar);
  HechtVelocityBoundary3d<CollisionOperator,Lattice,DevFlowYZ> vb(vel);
  // HechtVelocityBoundary3d<Lattice,VelocityOp> vb(vel);
  ierr = solver.addBoundaryCondition(Box3d(0,0,1,ny-2,1,nz-2),vb); CHKERRQ(ierr);

  PressureOp pr;
  pr.rho0 = 1.;
  HechtPressureBoundary3d<CollisionOperator,Lattice,PressureOp> pb(pr);
  // ierr = solver.addBoundaryCondition(Box3d(nx-1,nx-1,th+1,2*th,1,nz-2),pb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box3d(nx-1,nx-1,th+1,2*th,th+1,2*th),pb); CHKERRQ(ierr);

  /* Initial conditions */
  IsothermalMacros3d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  mac.uz = 0.;
  ierr = solver.uniformMacroInitialization(mac); CHKERRQ(ierr);

  IsothermalMacros3d*** macArray;
  Box3d arrayDelimiter;
  Box3d inflow(0,0,1,ny-2,1,nz-2);
  PetscScalar ux,uy,uz;

  ierr = solver.getInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);

  for (auto kk : arrayDelimiter.zRange){
    for (auto jj : arrayDelimiter.yRange){
      for (auto ii : arrayDelimiter.xRange){
        vel(0,ii,jj,kk,ux,uy,uz);
        macArray[kk][jj][ii].ux = ux;
        macArray[kk][jj][ii].uy = uy;
        macArray[kk][jj][ii].uz = uz;
      }
    }
  }

  ierr = solver.restoreInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);

  ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
  ierr = solver.collideAndSwap(); CHKERRQ(ierr);

  PetscInt numSteps = 1000;
  PetscInt step = 0;
  char output[50];
  double t1,t2;
  t1 = MPI_Wtime();
  while (step < numSteps){
    ierr = solver.streamBySwapping(); CHKERRQ(ierr);
    if (step % 50 == 0){
      sprintf(output,"Macros_step_%d.vts",step);
      ierr = solver.outputMacros(outputFolder,output); CHKERRQ(ierr);
    }
    // if (step % 5 == 0 && step <= 200){
    //   sprintf(output,"Dist_step_%d.vts",step);
    //   ierr = solver.outputDistributions(outputFolder,output); CHKERRQ(ierr);
    // }
    ++step;
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);
  }
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD,"Time: %f\n",t2-t1);
}
