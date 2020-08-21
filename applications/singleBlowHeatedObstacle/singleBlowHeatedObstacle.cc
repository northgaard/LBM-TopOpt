#include "TOpTEnLB2d.hh"
#include <memory>

static char help [] = "Implements a single blow optimization problem with heating solid.\n";

template <class Lattice, class Real>
// using IsoCollision = IncompressibleBGK2d<Lattice,Real>;
using IsoCollision = Cascaded2d<Lattice,Real>;
template <class Lattice, class Real>
using ThermalCollision = ThermalBGK2d<Lattice,Real>;
template <class Lattice, class Input>
// using VelocityBoundary = IncompressibleZouHeVelocity2d<Lattice,Input>;
using VelocityBoundary = StandardZouHeVelocity2d<Lattice,Input>;
template <class Lattice>
// using PressureBoundary = IncompressibleZouHePressure2d<Lattice>;
using PressureBoundary = StandardZouHePressure2d<Lattice>;

using IsoLattice = D2Q9;
using ThermalLattice = TD2Q4;

template <class Real>
using IsoInterpolation = BPInterpolation<Real>;
template <class Real>
using ThermalInterpolation = RampInterpolation<Real>;

struct RampUp {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy,
		  PetscScalar& ux, PetscScalar& uy) const
  {

    PetscScalar mult;

    if (t < rampUpTime){
      mult = sin(M_PI*t/(2.*rampUpTime));
    } else if (t < (rampUpTime + constTime)){
      mult = 1.;
    } else if (t < (rampUpTime + constTime + rampDownTime)){
      mult = cos(M_PI*(t - constTime - rampUpTime)/(2.*rampDownTime));
    } else {
      mult = 0.;
    }

    ux = mult*(vmax*(mid*mid - (idy - mid)*(idy - mid)))/(mid*mid);
    uy = 0.;

  }

  PetscScalar vmax;
  PetscScalar mid;
  PetscInt rampUpTime;
  PetscInt constTime;
  PetscInt rampDownTime;

};

struct NoRampUp {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy,
		  PetscScalar& ux, PetscScalar& uy) const
  {
    
    PetscScalar mult;
    if (t < constTime){
      mult = 1.;
    } else if (t < (constTime + rampTime)) {
      mult = cos(M_PI*(t - constTime)/(2.*rampTime));
    } else {
      mult = 0.;
    }
    ux = mult*(vmax*(mid*mid - (idy - mid)*(idy - mid)))/(mid*mid);
    uy = 0.;
    
  }

  PetscScalar vmax;
  PetscScalar mid;
  PetscInt constTime;
  PetscInt rampTime;
  
};

using InflowVelocity = NoRampUp;
// using InflowVelocity = RampUp;

struct ConstTemp {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy,
		  PetscScalar& T) const
  {
    T = T0;
  }

  PetscScalar T0;

};

PetscErrorCode setUpSolvers(const ObstacleLBSolverInfo2d& info,
			    PetscInt inflowLength,
			    ObstacleLBSolver2d& solver,
			    AdjointLBSolver2d& adjSolver)
{

  using Operator = HeatingThermalPartialBounceback2d
    <IsoLattice,ThermalLattice,IsoInterpolation,
     ThermalInterpolation,IsoCollision,ThermalCollision>;

  PetscErrorCode ierr;
  PetscInt nx = info.nx, ny = info.ny;
  PetscInt th = ny/3;

  ierr = solver.initializeSolver(info); CHKERRQ(ierr);
  ierr = adjSolver.initializeFromForwardSolver(solver); CHKERRQ(ierr);

  IncompressibleFlowParameters incPar;
  // incPar.lengthChar = th;
  incPar.lengthChar = ny-2;
  incPar.velocityChar = 0.1;
  incPar.ReynoldsNumber = 5000;

  ThermalFlowParameters thermPar;
  thermPar.incPar = incPar;
  // thermPar.PrandtlNumber = 0.71;
  thermPar.PrandtlNumber = 7.;

  // Set up interpolation
  PetscScalar isoPenal = 1.;
  IsoInterpolation<PetscScalar> isoInterp(isoPenal);

  PetscScalar Ck = 1e-2, thermPenal = 1.;
  ThermalInterpolation<PetscScalar> thermInterp(thermPenal,Ck);

  PetscScalar heatingCoefficient = 1e-2;
  // PetscScalar heatingCoefficient = 1.;
  Operator op(thermPar,heatingCoefficient,isoInterp,thermInterp);

  ObstacleCollisionLoop2d<Operator> colLoop(solver,op);
  solver.setCollisionFunction(colLoop);
  adjSolver.setAdjointCollisionFunction(colLoop.getAdjoint(solver));

  ThermalStreamingLoop2d<IsoLattice,ThermalLattice>
    strLoop(solver);
  solver.setStreamingFunction(strLoop);
  adjSolver.setAdjointStreamingFunction(strLoop.getAdjoint(solver));

  using IsoEq = typename IsoCollision<IsoLattice,PetscScalar>::
    Equilibrium;
  using ThermalEq = typename ThermalCollision<ThermalLattice,PetscScalar>::
    Equilibrium;

  ThermalEquilibriumInitializationLoop2d<IsoLattice,IsoEq,ThermalEq>
    eqInit(solver.getLocalBoundingBox());

  solver.setInitializationFunction(eqInit);

  ThermalMacros2d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  mac.T = 0.;
  solver.uniformMacroInitialization(mac);

  InflowVelocity vf;
  vf.vmax = incPar.velocityChar;
  vf.mid = (PetscScalar) (ny-1)/2.;
  // vf.rampUpTime = 100;
  // vf.constTime = 100;
  // vf.rampDownTime = 50;
  vf.constTime = 100;
  vf.rampTime = 100;

  ThermalMacros2d** macArray;
  Box2d inflow(0,0,1,ny-2);
  Box2d arrayDelimiter;
  PetscScalar ux,uy;

  solver.getInitialMacroArray(inflow,arrayDelimiter,&macArray);

  for (auto jj : arrayDelimiter.yRange){
    for (auto ii : arrayDelimiter.xRange){
      vf(0,ii,jj,ux,uy);
      macArray[jj][ii].ux = ux;
      macArray[jj][ii].uy = uy;
    }
  }

  solver.restoreInitialMacroArray(inflow,arrayDelimiter,&macArray);

  // Box2d initDomain(0,inflowLength,0,ny-1);

  // solver.getInitialMacroArray(initDomain,arrayDelimiter,&macArray);

  // for (auto jj : arrayDelimiter.yRange){
  //   for (auto ii : arrayDelimiter.xRange){
  //     macArray[jj][ii].T = -1.0;
  //   }
  // }

  // solver.restoreInitialMacroArray(inflow,arrayDelimiter,&macArray);

  auto bb = OnGridBounceback2d<IsoLattice>::make();
  addBoundaryCondition(Box2d(0,nx-1,0,0),bb,solver,adjSolver);
  addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),bb,solver,adjSolver);

  auto zhv = VelocityBoundary<IsoLattice,InflowVelocity>::
    make(vf);
  addBoundaryCondition(Box2d(0,0,1,ny-2),zhv,solver,adjSolver);

  PetscScalar rho0 = 1.;
  auto zhp = PressureBoundary<IsoLattice>::make(rho0);
  addBoundaryCondition(Box2d(nx-1,nx-1,1,ny-2),zhp,solver,adjSolver);

  auto tbb = ThermalBounceback2d<IsoLattice,ThermalLattice>::make();
  addBoundaryCondition(Box2d(0,nx-1,0,0),tbb,solver,adjSolver);
  addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),tbb,solver,adjSolver);
  // addBoundaryCondition(Box2d(nx-1,nx-1,1,ny-2),tbb,solver,adjSolver);
  auto tneu = ThermalZouHeNeumann2d<IsoLattice,ThermalLattice>::make();
  addBoundaryCondition(Box2d(nx-1,nx-1,1,ny-2),tneu,solver,adjSolver);
  // addBoundaryCondition(Box2d(0,0,1,ny-2),tneu,solver,adjSolver);

  ConstTemp ct;
  ct.T0 = 0.;
  auto tzh = ThermalZouHe2d<IsoLattice,ThermalLattice,ConstTemp>::
    make(ct);
  addBoundaryCondition(Box2d(0,0,1,ny-2),tzh,solver,adjSolver);

  // Box2d obst(inflowLength,inflowLength+th-1,th,2*th-1);
  // solver.addGeometricObstacle(obst,0.);

  return 0;

}

PetscErrorCode forwardSimulation(ObstacleLBSolver2d&,PetscInt);
PetscErrorCode optimizationLoop(UnsteadyLBTopOpt2d*,MMA*);

int main(int argc, char *argv[])
{

  TopTenInit init(argc,argv,help);

  PetscErrorCode ierr;

  ObstacleLBSolverInfo2d info;
  ObstacleLBSolver2d solver;
  AdjointLBSolver2d adjSolver;

  info.nx = 300;
  info.ny = 102;
  info.latticeDOF = IsoLattice::numDOF + ThermalLattice::numDOF;
  info.macroDOF = 4;
  PetscInt inflowLength = 80;
  PetscInt designLength = 30;
  // PetscInt th = info.ny/3;

  ierr = setUpSolvers(info,inflowLength,solver,adjSolver); CHKERRQ(ierr);

  PetscInt nt = 500;

  Box2d designDomain(inflowLength,inflowLength+designLength-1,1,info.ny-2);
  Box2d targetDomain(inflowLength+designLength,
		     inflowLength+designLength+20,1,info.ny-2);
  auto objective = DomainTemperature2d<IsoLattice,ThermalLattice>::
    make(targetDomain,solver,nt);

  // UnsteadyLBTopOpt2d top(solver,adjSolver,objective,nt,
  // 			 DensityFilter2d::make(4,designDomain,solver));
  UnsteadyLBTopOpt2d top(solver,adjSolver,objective,nt);
  top.setDesignDomain(designDomain);
  top.addVolumeConstraint(0.2,ConstrainMaterial);

  top.allocateMemoryForForwardSolution();

  auto mma = top.getMMA();

  // ierr = optimizationLoop(&top,mma.get()); CHKERRQ(ierr);

  // ierr = forwardSimulation(solver,nt); CHKERRQ(ierr);

  // double t1,t2;

  // t1 = MPI_Wtime();

  // top.computeObjective();

  // t2 = MPI_Wtime();

  // PetscPrintf(PETSC_COMM_WORLD,"Objective: %f\n",top.getObjective());
  // PetscPrintf(PETSC_COMM_WORLD,"Time for forward: %f\n",t2-t1);

  FiniteDifferenceCheck2d check(top,1e-4);

  // check.checkDomain(designDomain);
  Box2d dom(0,10,0,info.ny-1);
  check.checkDomain(dom);
  // DensityFilter2d filter(8,designDomain,solver);
  // DensityFilter2d filter(8,solver.getBoundingBox(),solver);

  // solver.outputDomain("Domain.vts");

  // filter.filterDesign(solver.obstacleGlobal,top.physicalFullDomainVec);
  // solver.setDomain(top.physicalFullDomainVec);

  // solver.outputDomain("FilteredDomain.vts");
  
  top.deallocateMemoryForForwardSolution();

}

PetscErrorCode optimizationLoop(UnsteadyLBTopOpt2d* top,
				MMA* mma)
{

  PetscErrorCode ierr;
  double t1,t2;
  PetscScalar ch = 1.;
  PetscScalar tol = 1e-4;
  PetscInt itr = 0;
  PetscInt maxItr = 1;
  Vec design;
  Vec designSens;
  PetscScalar obj,con;
  std::vector<Vec> consSens;
  PetscScalar movlim = 0.2;

  while (itr < maxItr && ch > tol){

    if ((itr % 10) == 0){
      top->outputFullDomain(itr);
    }

    t1 = MPI_Wtime();

    ierr = top->forwardAndAdjoint(); CHKERRQ(ierr);

    t2 = MPI_Wtime();

    ierr = top->getObjectiveConstraintsAndSensitivities(obj,con,
    						 designSens,
    						 consSens);
    CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"It.: %i, obj: %f, volume: %f, time for analysis: %f\n",
    		itr,obj,con,t2-t1);

    top->getDesignVector(design);

    mma->Update(design,designSens,&con,consSens.data(),movlim);
    mma->DesignChange(design,ch);

    PetscPrintf(PETSC_COMM_WORLD,"Design change: %f\n",ch);

    top->restoreDesignVector(design);

    ++itr;

  }

  top->outputFullDomain(itr);

  return 0;
  
}

PetscErrorCode forwardSimulation(ObstacleLBSolver2d& solver,
				 PetscInt numTimesteps)
{

  double t1,t2;

  PetscInt timestep = 0;

  solver.initializeDistributions();
  solver.outputMacros();
  // solver.outputDomain("Domain.vts");

  t1 = MPI_Wtime();

  while (timestep < numTimesteps){

    solver.streamAndCollide();
    ++timestep;

    if (timestep % 10 == 0){
      solver.outputMacros();
      // solver.outputDistributions();
    }
  }

  t2 = MPI_Wtime();

  PetscPrintf(PETSC_COMM_WORLD,"Time: %f\n",t2-t1);

  return 0;

}
