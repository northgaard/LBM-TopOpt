#include "TOpTEnLB2d.hh"

static char help [] = "Implementing the pump TopOpt problem.\n";
static char outputFolder[] = "/home/sebnorg/WorkData/NewCodePump/";

using Lattice = D2Q9;
// using BaseCollision = IncompressibleBGK2d<Lattice>;
using BaseCollision = IncompressibleMRTForcing2d<Lattice>;
using Interpolation = BPInterpolation<>;
using CollisionOperator = PartialBouncebackCollision<BaseCollision,Interpolation>;

// Boundary velocity
struct InflowVelocity {

  void operator()(PetscInt timestep, PetscInt idx, PetscInt idy,
                  PetscScalar& ux, PetscScalar& uy) const
  {
    ux = 0.;
    uy = -sin(M_PI*timestep/freq)*(vmax*(inletRad*inletRad - (idx - mid)*(idx - mid)))
      / (inletRad*inletRad);
  }
  PetscScalar vmax;
  PetscScalar mid;
  PetscScalar inletRad;
  PetscInt freq;
};

struct ConstPres {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy, PetscScalar& rho) const
  {
    rho = rho0;
  }
  PetscScalar rho0;
};

int main(int argc, char *argv[])
{
  TopTenInit topinit(argc,argv,help);
  PetscErrorCode ierr;
  IsothermalLBSolver2d<CollisionOperator> solver;

  NewLBSolverInfo2d info;
  PetscInt nx = 234;
  PetscInt ny = 234;
  PetscInt channelLength = 55;
  info.nx = nx;
  info.ny = ny;
  PetscInt th = ny/3;

  IncompressibleFlowParameters par;
  par.velocityChar = 0.05;
  par.ReynoldsNumber = 80.;
  par.lengthChar = th;
  PetscScalar nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
  PetscScalar ug = 0.5*par.velocityChar;
  PetscScalar inletDiam = th;
  PetscScalar g = 8.*nu*ug/(inletDiam*inletDiam);

  BaseCollision baseOp(par,-g,0.);
  Interpolation interp(0.5);
  CollisionOperator colOp(baseOp,interp);
  solver.make(PETSC_COMM_WORLD,info,colOp);

  /* Define boundaries */
  OnGridBounceback2d<Lattice> bb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,0,1,th),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,1,th),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,0,2*th+1,ny-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,2*th+1,ny-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(1,th,ny-1,ny-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(2*th+1,nx-2,ny-1,ny-1),bb); CHKERRQ(ierr);

  PetscScalar rho0 = 1.;
  ConstPres outPres;
  outPres.rho0 = rho0;
  // PetscScalar delta = 0.;
  // ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhp(outPres);
  // ierr = solver.addBoundaryCondition(Box2d(0,0,th+1,2*th),zhp); CHKERRQ(ierr);
  // ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,th+1,2*th),zhp); CHKERRQ(ierr);
  ZouHeNeumann2d<CollisionOperator,Lattice> zhN;
  ierr = solver.addBoundaryCondition(Box2d(0,0,th+1,2*th),zhN); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,th+1,2*th),zhN); CHKERRQ(ierr);

  InflowVelocity vf;
  // vf.vmax = par.velocityChar;
  vf.vmax = 0.;
  vf.mid = (PetscScalar) (ny-1)/2.;
  vf.inletRad = (PetscScalar) th/2.;
  vf.freq = 400;
  ZouHeVelocity2d<CollisionOperator,Lattice,InflowVelocity> zhv(vf);
  ierr = solver.addBoundaryCondition(Box2d(th+1,2*th,ny-1,ny-1),zhv); CHKERRQ(ierr);

  /* Initial conditions */
  IsothermalMacros2d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  ierr = solver.uniformMacroInitialization(mac); CHKERRQ(ierr);

  /* Place material */
  PetscInt lowerHeight = 10;
  PetscScalar init[1] = {0.};
  ierr = solver.setUniformFieldValues(Box2d(0,nx-1,0,lowerHeight-1),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(0,channelLength-1,0,th),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(nx-channelLength,nx-1,0,th),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(0,channelLength-1,2*th+1,ny-1),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(nx-channelLength,nx-1,2*th+1,ny-1),init);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(0,th,ny-channelLength,ny-1),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(2*th+1,nx-1,ny-channelLength,ny-1),init);
  CHKERRQ(ierr);

  ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
  ierr = solver.collideAndSwap(); CHKERRQ(ierr);

  AdjointLBSolver adjSolver;
  ierr = solver.makeSourceAdjoint(&adjSolver); CHKERRQ(ierr);


  Vec fullDomain, filteredDomain, designDomain;
  FullToDesignMapping2d map(Box2d(channelLength,nx-1-channelLength,
                                  lowerHeight,ny-1-channelLength),solver);
  ProjectionFilter2d filter(5,Box2d(channelLength,nx-1-channelLength,
                                    lowerHeight,ny-1-channelLength),
                            solver,adjSolver);
  filter.setBeta(1.);
  filter.setEta(0.1);
  ierr = solver.createDesignVec(&fullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&filteredDomain); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designDomain); CHKERRQ(ierr);
  ierr = VecSet(designDomain,1.); CHKERRQ(ierr);
  ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
  ierr = filter.filterDesign(fullDomain,filteredDomain); CHKERRQ(ierr);
  ierr = solver.setFieldFromDesignVec(filteredDomain); CHKERRQ(ierr);

  ierr = solver.outputFields(outputFolder,"PumpLayout.vts"); CHKERRQ(ierr);

  PetscInt numSteps = 24000;
  PetscInt step = 0;
  double t1,t2;
  OutflowEast2d<CollisionOperator> objective(Box2d(nx-1,nx-1,th+1,2*th),
                                             solver);
  // NaiveCheckpointing check(numSteps+1,solver.distributionsLocal);
  LBStaticCheckpointing check(500,numSteps+1,solver);
  PetscScalar timeAveraging = (PetscScalar) 1./numSteps;
  PetscScalar totalOut = 0.,tOut;
  char output[50];

  t1 = MPI_Wtime();

  while (step < numSteps){
    ierr = solver.streamBySwapping(); CHKERRQ(ierr);
    ++step;
    ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
    solver.collideAndSwap(); CHKERRQ(ierr);
    ierr = solver.computeMacroFunctional(objective,timeAveraging,&tOut); CHKERRQ(ierr);
    totalOut += tOut;

    if (step % 500 == 0){
      sprintf(output,"Macros_step_%d.vts",step);
      ierr = solver.outputMacros(outputFolder,output); CHKERRQ(ierr);
    }
  }

  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD,"Outflow: %f\n",totalOut);
  PetscPrintf(PETSC_COMM_WORLD,"Forward time: %f\n",t2-t1);

  ierr = VecDestroy(&fullDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&filteredDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&designDomain); CHKERRQ(ierr);
}
