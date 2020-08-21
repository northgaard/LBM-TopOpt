#include "TOpTEnLB2d.hh"

static char help [] = "Pressure driven diode for the AD paper.\n";
static char outputFolder[] = "/home/sebnorg/WorkData/PressureDiode/";

using Lattice = D2Q9;
using BaseCollision = IncompressibleMRT2d<Lattice>;
// using BaseCollision = IncompressibleBGK2d<Lattice>;
// using BaseCollision = Cascaded2d<Lattice>;
using Interpolation = BPInterpolation<>;
using CollisionOperator = PartialBouncebackCollision<BaseCollision,Interpolation>;

struct ConstPres {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy, PetscScalar& rho) const
  {
    rho = rho0;
  }
  PetscScalar rho0;
};

struct OscPres {

  void operator()(PetscInt t, PetscInt, PetscInt, PetscScalar& rho) const
  {
    PetscInt denom = fullPeriod/2;
    rho = rho0 + deltaRho*sin(M_PI*t/denom);
  }
  PetscScalar rho0;
  PetscScalar deltaRho;
  PetscInt fullPeriod;
};

struct InflowPres {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy, PetscScalar& rho) const
  {
    if (t < upTime){
      rho = rho0 + deltaRho;
    } else if (upTime <= t && t < upTime+linTime){
      rho = rho0 + deltaRho - (t-upTime)*(deltaRho/linTime);
    } else if (t >= upTime+linTime){
      rho = rho0;
    }
  }
  PetscScalar rho0;
  PetscScalar deltaRho;
  PetscInt upTime,linTime;
};

struct OutflowPres {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy, PetscScalar& rho) const
  {
    if (t < upTime){
      rho = rho0;
    } else if (upTime <= t && t < upTime+linTime){
      rho = rho0 + (t-upTime)*(deltaRho/linTime);
    } else if (t >= upTime+linTime){
      rho = rho0 + deltaRho;
    }
  }
  PetscScalar rho0;
  PetscScalar deltaRho;
  PetscInt upTime,linTime;
};

int main(int argc, char *argv[])
{
  TopTenInit topinit(argc,argv,help);
  PetscErrorCode ierr;
  IsothermalLBSolver2d<CollisionOperator> solver;

  NewLBSolverInfo2d info;
  PetscInt nx = 350;
  PetscInt ny = 125;
  PetscInt inflowLength = 25;
  info.nx = nx;
  info.ny = ny;
  PetscInt fi = ny/5;

  IncompressibleFlowParameters par;
  par.velocityChar = 0.1;
  par.ReynoldsNumber = 10;
  par.lengthChar = fi;

  BaseCollision baseOp(par);
  Interpolation interp(1.);
  CollisionOperator colOp(baseOp,interp);
  solver.make(info,colOp);

  /* Define boundaries */
  OnGridBounceback2d<Lattice> bb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,0,1,2*fi),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,0,3*fi+1,ny-2),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,1,2*fi),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,3*fi+1,ny-2),bb); CHKERRQ(ierr);

  // PetscScalar rho0 = 1.;
  // PetscScalar deltaRho = 0.1;
  PetscInt totSteps = 20000;
  // PetscInt linTime = 1000;
  // ConstPres inPres;
  // inPres.rho0 = rho0 + deltaRho;
  // ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhpIn(inPres);
  // InflowPres inPres;
  // inPres.rho0 = rho0;
  // inPres.deltaRho = deltaRho;
  // inPres.upTime = (totSteps - linTime)/2;
  // inPres.linTime = linTime;
  // ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,InflowPres> zhpIn(inPres);
  // ierr = solver.addBoundaryCondition(Box2d(0,0,2*fi+1,3*fi),zhpIn); CHKERRQ(ierr);
  // ConstPres outPres;
  // outPres.rho0 = rho0;
  // ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhpOut(outPres);
  // OutflowPres outPres;
  // outPres.rho0 = rho0;
  // outPres.deltaRho = deltaRho;
  // outPres.upTime = inPres.upTime;
  // outPres.linTime = linTime;
  // ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,OutflowPres> zhpOut(outPres);
  // ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,2*fi+1,3*fi),zhpOut); CHKERRQ(ierr);
  OscPres inPres;
  inPres.rho0 = 1.;
  inPres.deltaRho = 0.1;
  inPres.fullPeriod = 4000;
  ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,OscPres> zhpIn(inPres);
  ierr = solver.addBoundaryCondition(Box2d(0,0,2*fi+1,3*fi),zhpIn); CHKERRQ(ierr);
  ConstPres outPres;
  outPres.rho0 = 1.;
  ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhpOut(outPres);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,2*fi+1,3*fi),zhpOut); CHKERRQ(ierr);

  /* Initial conditions */
  IsothermalMacros2d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  ierr = solver.uniformMacroInitialization(mac); CHKERRQ(ierr);
  // IsothermalMacros2d** macArray;
  // Box2d inflow(0,0,2*fi+1,3*fi);
  // Box2d arrayDelimiter;

  // ierr = solver.getInitialMacroArray(inflow,&arrayDelimiter,
  //                                    &macArray); CHKERRQ(ierr);

  // for (auto jj : arrayDelimiter.yRange){
  //   for (auto ii : arrayDelimiter.xRange){
  //     macArray[jj][ii].rho = rho0 + deltaRho;
  //   }
  // }

  // ierr = solver.restoreInitialMacroArray(inflow,&arrayDelimiter,
  //                                         &macArray); CHKERRQ(ierr);

  /* Place material */
  PetscScalar init[1] = {0.};
  ierr = solver.setUniformFieldValues(Box2d(0,inflowLength-1,1,2*fi),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(0,inflowLength-1,3*fi+1,ny-2),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(nx-inflowLength,nx-1,1,2*fi),init); CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(nx-inflowLength,nx-1,3*fi+1,ny-2),init);
  CHKERRQ(ierr);

  AdjointLBSolver adjSolver;
  ierr = solver.makeSourceAdjoint(&adjSolver); CHKERRQ(ierr);

  Vec fullDomain, designDomain;
  Vec physicalFullDomain, physicalDesignDomain;
  Vec sens,designSens;
  Vec vcSens,designVcSens;
  FullToDesignMapping2d map(Box2d(inflowLength,nx-inflowLength-1,1,ny-2),solver);
  DensityFilter2d filter(3,Box2d(inflowLength,nx-inflowLength-1,1,ny-2),
                         solver,adjSolver);

  ierr = solver.createDesignVec(&fullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&physicalFullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&vcSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designDomain); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designVcSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&physicalDesignDomain); CHKERRQ(ierr);
  ierr = VecSet(designDomain,1.); CHKERRQ(ierr);

  // ierr = solver.outputFields(outputFolder,"Init.vts"); CHKERRQ(ierr);

  char output[50];

  OutflowEast2d<CollisionOperator> objective(Box2d(nx-1,nx-1,2*fi+1,3*fi),
                                             solver.getLocalBoundingBox());
  // PositiveUxFlow2d<CollisionOperator> objective(Box2d(0,0,2*fi+1,3*fi),
  //                                               Box2d(nx-1,nx-1,2*fi+1,3*fi),
  //                                               solver.getLocalBoundingBox());
  VolumeConstraint vc(0.5,ConstrainFluid,designDomain);
  PetscScalar volumeValue;
  // ierr = map.mapFullToDesign(physicalFullDomain,physicalDesignDomain); CHKERRQ(ierr);
  // vc.computeFunctional(physicalDesignDomain,&volumeValue);
  // PetscPrintf(PETSC_COMM_WORLD,"Volume constraint: %f\n",volumeValue);
  LBStaticCheckpointing check(1000,totSteps+1,solver);
  // NaiveCheckpointing check(totSteps+1,solver.distributionsLocal);
  PetscScalar timeAveraging = (PetscScalar) 1./totSteps;
  PetscScalar tstepOut;
  PetscScalar finalObj;
  PetscInt numConstraints = 1;
  PetscInt numDesign;
  ierr = VecGetSize(designDomain,&numDesign); CHKERRQ(ierr);
  MMA mma(numDesign,numConstraints,designDomain);
  PetscScalar tol = 1e-4;
  PetscScalar ch = 1e-5;
  PetscInt maxIter = 1000;
  PetscInt curIter = 0;
  PetscScalar movlim = 0.2;
  PetscInt step;
  PetscPrintf(PETSC_COMM_WORLD,"-----------------------\n");

  while (ch > tol && curIter < maxIter){

    ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
    ierr = filter.filterDesign(fullDomain,physicalFullDomain); CHKERRQ(ierr);
    ierr = solver.setFieldFromDesignVec(physicalFullDomain); CHKERRQ(ierr);
    check.reset(nullptr);

    if (curIter % 50 == 0){
      sprintf(output,"Design_iter_%d.vts",curIter);
      ierr = solver.outputFields(outputFolder,output); CHKERRQ(ierr);
    }

    ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
    ierr = check.makeCheckpoint(0,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);

    double t1,t2;

    t1 = MPI_Wtime();
    step = 0;
    finalObj = 0.;

    while (step < totSteps){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      // if (step % 1000 == 0){
      //   sprintf(output,"Macros_step_%d.vts",step);
      //   ierr = solver.outputMacros(outputFolder,output); CHKERRQ(ierr);
      // }
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
      ierr = solver.computeMacroFunctional(objective,timeAveraging,&tstepOut); CHKERRQ(ierr);
      finalObj += tstepOut;
    }

    t2 = MPI_Wtime();

    PetscPrintf(PETSC_COMM_WORLD,"Forward time: %f\n",t2-t1);

    ierr = map.mapFullToDesign(physicalFullDomain,physicalDesignDomain); CHKERRQ(ierr);
    vc.computeFunctional(physicalDesignDomain,&volumeValue);
    vc.computeSensitivities(vcSens);

    ierr = adjSolver.resetAdjoints(); CHKERRQ(ierr);
    ierr = adjSolver.setCurrentTimestep(step);

    t1 = MPI_Wtime();

    while (step > 0){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal,
                                      objective,timeAveraging); CHKERRQ(ierr);
      ierr = adjSolver.adjointStream(objective); CHKERRQ(ierr);
      --step;
    }
    ierr = adjSolver.getSensitivities(&sens); CHKERRQ(ierr);
    // VecView(sens,PETSC_VIEWER_STDOUT_WORLD);

    t2 = MPI_Wtime();

    PetscPrintf(PETSC_COMM_WORLD,"Adjoint time: %f\n",t2-t1);

    ierr = filter.filterSensitivities(sens,&vcSens,numConstraints); CHKERRQ(ierr);
    // Map to design domain
    ierr = map.mapFullToDesign(sens,designSens); CHKERRQ(ierr);
    ierr = map.mapFullToDesign(vcSens,designVcSens); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"Iter: %d, objective: %f, volume constraint: %f\n",
                curIter,finalObj,volumeValue);

    mma.Update(designDomain,designSens,&volumeValue,&designVcSens,movlim);
    mma.DesignChange(designDomain,ch);
    ++curIter;
    PetscPrintf(PETSC_COMM_WORLD,"Design change: %f\n",ch);
    PetscPrintf(PETSC_COMM_WORLD,"-----------------------\n");
  }

  // Output final design
  ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
  ierr = filter.filterDesign(fullDomain,physicalFullDomain); CHKERRQ(ierr);
  ierr = solver.setFieldFromDesignVec(physicalFullDomain); CHKERRQ(ierr);
  ierr = solver.outputFields(outputFolder,"Final_design.vts"); CHKERRQ(ierr);
}
