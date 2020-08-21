#include "TOpTEnLB2d.hh"

static char help [] = "Implementing the robust pump TopOpt problem.\n";
static char outputFolder[] = "/home/sebnorg/";

using Lattice = D2Q9;
// using BaseCollision = IncompressibleBGK2d<Lattice>;
using BaseCollision = IncompressibleMRT2d<Lattice>;
using Interpolation = BPInterpolation<>;
using CollisionOperator = PartialBouncebackCollision<BaseCollision,Interpolation>;

// Boundary velocity
struct InflowVelocity {
  void operator()(PetscInt timestep, PetscInt idx, PetscInt idy,
                  PetscScalar& ux, PetscScalar& uy) const
  {
    ux = 0.;
    uy = -sin(2.*M_PI*timestep/freq)*(vmax*(inletRad*inletRad - (idx - mid)*(idx - mid)))
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
  constexpr PetscMPIInt numCases = 3;
  MPIEvenSplit split(numCases);
  MPI_Comm splitComm;
  PetscMPIInt splitID;
  PetscMPIInt* roots;
  split.getCommunicationInfo(&splitComm,&splitID,&roots);
  // PetscScalar etaValues[numCases] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  PetscScalar etaValues[numCases] = {0.3,0.5,0.7};

  NewLBSolverInfo2d info;
  PetscInt nx = 234;
  PetscInt ny = 234;
  PetscInt channelLength = 50;
  info.nx = nx;
  info.ny = ny;
  PetscInt th = ny/3;

  PetscScalar dx = 1./(th-1);
  PetscScalar dt = dx*dx;
  PetscInt cycle = std::round(1./dt);
  PetscScalar Re = 80.;
  PetscScalar nu = (dt/(dx*dx))*(1./Re);
  PetscScalar omega = 1./(Lattice::csSqInv*nu + 0.5);
  PetscScalar uc = dt/dx;
  PetscInt numSteps = 3*cycle;

  PetscPrintf(PETSC_COMM_WORLD,"dx: %f, dt: %f\n",dx,dt);
  PetscPrintf(PETSC_COMM_WORLD,"Length of cycle: %d, number of time steps: %d\n",cycle,numSteps);
  PetscPrintf(PETSC_COMM_WORLD,"Reynolds number: %f, omega: %f\n",Re,omega);
  PetscPrintf(PETSC_COMM_WORLD,"Maximum inflow velocity (lu): %f\n",uc);

  BaseCollision baseOp(omega);
  Interpolation interp(0.5);
  CollisionOperator colOp(baseOp,interp);
  solver.make(splitComm,info,colOp);

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
  ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhp(outPres);
  ierr = solver.addBoundaryCondition(Box2d(0,0,th+1,2*th),zhp); CHKERRQ(ierr);
  PetscScalar dr = 0.;
  ConstPres incPres;
  incPres.rho0 = rho0 + dr;
  ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhpInc(incPres);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,th+1,2*th),zhpInc); CHKERRQ(ierr);

  InflowVelocity vf;
  vf.vmax = uc;
  vf.mid = (PetscScalar) (ny-1)/2.;
  vf.inletRad = (PetscScalar) th/2.;
  vf.freq = cycle;
  ZouHeVelocity2d<CollisionOperator,Lattice,InflowVelocity> zhv(vf);
  ierr = solver.addBoundaryCondition(Box2d(th+1,2*th,ny-1,ny-1),zhv); CHKERRQ(ierr);

  /* Initial conditions */
  IsothermalMacros2d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  ierr = solver.uniformMacroInitialization(mac); CHKERRQ(ierr);

  /* Place material */
  PetscInt lowerHeight = 20;
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

  Vec fullDomain, physicalDomain;
  Vec designDomain, designDomainWorld, physicalDesignDomain;
  Vec objDummySens;
  FullToDesignMapping2d map(Box2d(channelLength,nx-1-channelLength,
                                  lowerHeight,ny-1-channelLength),solver);
  // ProjectionFilter2d filter(12,Box2d(channelLength,nx-1-channelLength,
  //                                   lowerHeight,ny-1-channelLength),
  //                           solver,adjSolver);
  DoubleProjectionFilter2d filter(8,4,Box2d(channelLength,nx-1-channelLength,
                                    lowerHeight,ny-1-channelLength),
                            solver,adjSolver);
  filter.setFirstBeta(2.);
  filter.setSecondBeta(1.);
  filter.setFirstEta(0.5);
  filter.setSecondEta(etaValues[splitID]);
  ierr = solver.createDesignVec(&fullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&physicalDomain); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designDomain); CHKERRQ(ierr);
  ierr = map.createWorldDesignDomainVec(&designDomainWorld); CHKERRQ(ierr);
  ierr = map.createWorldDesignDomainVec(&objDummySens); CHKERRQ(ierr);
  ierr = VecSet(designDomainWorld,1.); CHKERRQ(ierr);
  ierr = VecSet(objDummySens,0.); CHKERRQ(ierr);

  ierr = split.setUpScatter(designDomain,designDomainWorld); CHKERRQ(ierr);

  char output[100];
  // sprintf(output,"PumpLayout_%d.vts",splitID);
  // ierr = solver.outputFields(outputFolder,output); CHKERRQ(ierr);

  PetscInt step = 0;
  double t1,t2,t3,t4;
  // OutflowEast2d<CollisionOperator> objective(Box2d(nx-1,nx-1,th+1,2*th),solver);
  ThroughFlow2d<CollisionOperator> objective(Box2d(0,0,th+1,2*th),Box2d(nx-1,nx-1,th+1,2*th),
                                             solver);
  LBStaticCheckpointing check(10,numSteps+1,solver);
  PetscScalar timeAveraging = (PetscScalar) 1./numSteps;
  PetscScalar totalOut,tOut;
  constexpr PetscInt numConstraints = numCases + 1;
  PetscScalar conValues[numConstraints];
  PetscInt numDesign;
  ierr = VecGetSize(designDomainWorld,&numDesign); CHKERRQ(ierr);
  MMA mma(numDesign,numConstraints,designDomainWorld);
  PetscScalar tol = 1e-4;
  PetscScalar ch = 1.;
  PetscInt maxIter = 1000;
  PetscScalar betaInc = 1.5;
  PetscInt incFreq = 50;
  PetscInt curIter = 0;
  PetscScalar movlim = 0.1;
  constexpr PetscScalar magicAdd = 5.;
  {
    PetscScalar robustConstants[numConstraints];
    for (size_t ii = 0; ii < numConstraints-1; ++ii){
      robustConstants[ii] = 1.;
    }
    robustConstants[numConstraints-1] = 0.;
    mma.SetRobustConstants(robustConstants);
  }

  /* Sensitivity vectors */
  Vec worldSens[numConstraints];
  Vec sens,designSens;
  Vec vcSens,designVcSens;

  for (PetscInt ii = 0; ii < numConstraints; ++ii){
    ierr = map.createWorldDesignDomainVec(&(worldSens[ii])); CHKERRQ(ierr);
  }
  ierr = map.createDesignDomainVec(&designSens); CHKERRQ(ierr);
  // ID of the collection of processors that will do volume constraint
  PetscInt volumeID = 0;
  if (splitID == volumeID){
    ierr = map.createDesignDomainVec(&designVcSens); CHKERRQ(ierr);
    ierr = map.createDesignDomainVec(&physicalDesignDomain); CHKERRQ(ierr);
    ierr = solver.createDesignVec(&vcSens); CHKERRQ(ierr);
  } else {
    vcSens = nullptr;
    designVcSens = nullptr;
    physicalDesignDomain = nullptr;
  }
  VolumeConstraint vc(0.8,ConstrainFluid,designDomain);
  PetscScalar volumeValue;
  PetscInt localConstraints = (splitID == volumeID ? 1 : 0);

  PetscPrintf(PETSC_COMM_WORLD,"--------------------------\n");

  while (ch > tol && curIter < maxIter){

    ierr = split.scatterCommonFromWorld(designDomainWorld,designDomain); CHKERRQ(ierr);
    ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
    ierr = filter.filterDesign(fullDomain,physicalDomain); CHKERRQ(ierr);
    ierr = solver.setFieldFromDesignVec(physicalDomain); CHKERRQ(ierr);
    check.reset(nullptr);

    if (curIter % 50 == 0){
      sprintf(output,"Design_iter_%d_realization_%d.vts",curIter,splitID);
      ierr = solver.outputFields(outputFolder,output); CHKERRQ(ierr);
      if (splitID == 0){
        sprintf(output,"Design_variables_iter_%d.bin",curIter);
        ierr = writeVectorToBinary(outputFolder,output,designDomain,splitComm);
      }
    }

    ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
    ierr = check.makeCheckpoint(0,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);

    t1 = MPI_Wtime();
    step = 0;
    totalOut = 0.;

    while (step < numSteps){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      solver.collideAndSwap(); CHKERRQ(ierr);
      ierr = solver.computeMacroFunctional(objective,timeAveraging,&tOut); CHKERRQ(ierr);
      totalOut += tOut;
    }

    t2 = MPI_Wtime();

    // Broadcast objective to world
    conValues[splitID] = totalOut;
    for (PetscInt ii = 0; ii < numCases; ++ii){
      MPI_Bcast(&(conValues[ii]),1,MPIU_SCALAR,roots[ii],PETSC_COMM_WORLD);
      PetscPrintf(PETSC_COMM_WORLD,"Objective, case %d: %f\n",ii,conValues[ii]);
      // Add the magic number!
      conValues[ii] += magicAdd;
    }

    /* Adjoint solve */
    ierr = adjSolver.resetAdjoints(); CHKERRQ(ierr);
    ierr = adjSolver.setCurrentTimestep(step); CHKERRQ(ierr);

    t3 = MPI_Wtime();

    while (step > 0){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,
                                      solver.materialLocal); CHKERRQ(ierr);
      ierr = adjSolver.computeCollideSource(objective,timeAveraging,solver.distributionsLocal,
                                            solver.materialLocal); CHKERRQ(ierr);
      ierr = adjSolver.adjointStream(solver.distributionsLocal); CHKERRQ(ierr);
      --step;
    }
    ierr = adjSolver.getSensitivities(&sens); CHKERRQ(ierr);

    t4 = MPI_Wtime();

    if (splitID == volumeID){
      ierr = map.mapFullToDesign(physicalDomain,physicalDesignDomain); CHKERRQ(ierr);
      vc.computeFunctional(physicalDesignDomain,&volumeValue);
      vc.computeSensitivities(vcSens);
      conValues[numCases] = volumeValue;
    }
    // Broadcast volume constraint value
    MPI_Bcast(&(conValues[numCases]),1,MPIU_SCALAR,roots[volumeID],PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD,"Volume constraint: %f\n",conValues[numCases]);
    PetscPrintf(PETSC_COMM_WORLD,"Forward time: %f\n",t2-t1);
    PetscPrintf(PETSC_COMM_WORLD,"Adjoint time: %f\n",t4-t3);

    ierr = filter.filterSensitivities(sens,&vcSens,localConstraints); CHKERRQ(ierr);
    ierr = map.mapFullToDesign(sens,designSens); CHKERRQ(ierr);

    if (splitID == volumeID){
      ierr = map.mapFullToDesign(vcSens,designVcSens); CHKERRQ(ierr);
    }

    // Scatter objective sensitivities to world
    ierr = split.scatterToWorld(designSens,worldSens); CHKERRQ(ierr);
    // Scatter volume constraint sensitivities to world
    ierr = split.scatterSingleToWorld(designVcSens,worldSens[numCases],volumeID);

    mma.Update(designDomainWorld,objDummySens,conValues,worldSens,movlim);
    mma.DesignChange(designDomainWorld,ch);
    ++curIter;
    if (curIter % incFreq == 0 && curIter <= 500){
      filter.increaseBeta(betaInc);
    }
    PetscPrintf(PETSC_COMM_WORLD,"Iteration: %d\n",curIter);
    PetscPrintf(PETSC_COMM_WORLD,"Design change: %f\n",ch);
    PetscPrintf(PETSC_COMM_WORLD,"--------------------------\n");
  }

  // Output final design
  ierr = split.scatterCommonFromWorld(designDomainWorld,designDomain); CHKERRQ(ierr);
  ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
  ierr = filter.filterDesign(fullDomain,physicalDomain); CHKERRQ(ierr);
  ierr = solver.setFieldFromDesignVec(physicalDomain); CHKERRQ(ierr);

  sprintf(output,"Final_design_realization_%d.vts",splitID);
  ierr = solver.outputFields(outputFolder,output); CHKERRQ(ierr);
  if (splitID == 0){
    ierr = writeVectorToBinary(outputFolder,"Design_variables_final.bin",designDomain,splitComm);
  }

  ierr = VecDestroy(&fullDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&physicalDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&designDomain); CHKERRQ(ierr);
}
