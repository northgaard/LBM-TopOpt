#include "TOpTEnLB2d.hh"

static char help [] = "Example: implementing the diffuser TopOpt problem.";
static char outputFolder[] = "/home/sebnorg/WorkData/DiffuserTest/";

using Lattice = D2Q9;
// using BaseCollision = IncompressibleBGK2d<Lattice>;
// using BaseCollision = IncompressibleMRT2d<Lattice>;
using BaseCollision = IncompressibleCascaded2d<Lattice>;
using Interpolation = BPInterpolation<>;
using CollisionOperator = PartialBouncebackCollision<BaseCollision,Interpolation>;

struct EtaComputation {
  PetscScalar operator()(PetscInt,PetscInt) const
  {
    return 0.5;
  }
};

// Boundary velocity class
struct InflowVelocity {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy,
                  PetscScalar& ux, PetscScalar& uy) const
  {
    ux = (vmax*(mid*mid - (idy - mid)*(idy - mid)))/(mid*mid);
    // ux = 0.;
    uy = 0.;
  }
  PetscScalar vmax;
  PetscScalar mid;
};

struct ConstPres {

  void operator()(PetscInt t, PetscInt idx, PetscInt idy,
                  PetscScalar& rho) const
  {
    rho = rho0;
  }
  PetscScalar rho0;
};

class ComputeObjective {

public:
  ComputeObjective(PetscInt _nt) : numSteps(_nt)
  {
    timeAveraging = (PetscScalar) 1./numSteps;
  }
  PetscErrorCode operator()(NewLBSolver2d& solver, const LBMacroFunctional& func,
                            PetscScalar* value, bool dumpMacros = false)
  {
    PetscErrorCode ierr;
    PetscInt step = 0;
    PetscScalar tvalue;
    PetscFunctionBeginUser;
    ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    *value = 0.;
    while (step < numSteps){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
      if (step % 50 == 0){
        char output[50];
        sprintf(output,"Macros_iter_%d.vts",step);
        ierr = solver.outputMacros(outputFolder,output); CHKERRQ(ierr);
      }
      ierr = solver.computeMacroFunctional(func,timeAveraging,&tvalue); CHKERRQ(ierr);
      *value += tvalue;
    }
    PetscFunctionReturn(0);
  }
private:
  PetscInt numSteps;
  PetscScalar timeAveraging;
};

class ComputeObjectiveAndSens {

public:
  ComputeObjectiveAndSens(CheckpointingBase& _chk, PetscInt _nt)
    : numSteps(_nt), check(_chk)
  {
    timeAveraging = (PetscScalar) 1./numSteps;
  }
  PetscErrorCode operator()(NewLBSolver2d& solver, AdjointLBSolver& adjSolver,
                            const LBMacroFunctional& func, PetscScalar* value,
                            Vec* sens)
  {
    PetscErrorCode ierr;
    PetscInt step = 0;
    PetscScalar tvalue;
    PetscFunctionBeginUser;
    check.reset(nullptr);
    ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
    ierr = check.makeCheckpoint(0,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    *value = 0.;
    while (step < numSteps){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
      ierr = solver.computeMacroFunctional(func,timeAveraging,&tvalue); CHKERRQ(ierr);
      *value += tvalue;
    }
    adjSolver.resetAdjoints();
    while (step > 0){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      // ierr = adjSolver.adjointCollide(solver.distributionsLocal,
      //                                 solver.materialLocal,
      //                                 func,timeAveraging); CHKERRQ(ierr);
      // ierr = adjSolver.adjointStream(func);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,
                                      solver.materialLocal); CHKERRQ(ierr);
      ierr = adjSolver.computeCollideSource(func,timeAveraging,solver.distributionsLocal,
                                            solver.materialLocal); CHKERRQ(ierr);
      ierr = adjSolver.adjointStream(solver.distributionsLocal); CHKERRQ(ierr);
      --step;
    }
    adjSolver.getSensitivities(sens);
    PetscFunctionReturn(0);
  }
private:
  PetscInt numSteps;
  CheckpointingBase& check;
  PetscScalar timeAveraging;
};

int main(int argc, char *argv[])
{
  TopTenInit init(argc,argv,help);
  PetscErrorCode ierr;
  IsothermalLBSolver2d<CollisionOperator> solver;

  NewLBSolverInfo2d info;
  PetscInt nx = 150;
  PetscInt ny = 102;
  info.nx = nx;
  info.ny = ny;
  PetscInt th = ny / 3;

  IncompressibleFlowParameters par;
  par.velocityChar = 0.1;
  par.ReynoldsNumber = 100.0;
  par.lengthChar = ny;

  BaseCollision baseOp(par);
  Interpolation interp(1.);
  CollisionOperator colOp(baseOp,interp);
  solver.make(PETSC_COMM_WORLD,info,colOp);

  /* Define boundaries */
  OnGridBounceback2d<Lattice> bb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,1,th),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,2*th+1,ny-2),bb); CHKERRQ(ierr);

  PetscScalar rho0 = 1.;
  ConstPres outP;
  outP.rho0 = rho0;
  ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhp(outP);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,th+1,2*th),zhp); CHKERRQ(ierr);
  // ZouHeNeumann2d<CollisionOperator,Lattice> zhN;
  // ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,th+1,2*th),zhN); CHKERRQ(ierr);

  // PetscScalar deltaRho = 0.1;
  // ConstPres inP;
  // inP.rho0 = rho0 + deltaRho;
  // ZouHePressure2d<Lattice,CollisionOperator::Equilibrium,ConstPres> zhpIn(inP);
  // ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),zhpIn); CHKERRQ(ierr);

  InflowVelocity vf;
  vf.vmax = par.velocityChar;
  vf.mid = (PetscScalar) (ny-1)/2.;
  ZouHeVelocity2d<CollisionOperator,Lattice,InflowVelocity> zhv(vf);
  ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),zhv); CHKERRQ(ierr);


  /* Initial conditions */
  IsothermalMacros2d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  ierr = solver.uniformMacroInitialization(mac); CHKERRQ(ierr);

  IsothermalMacros2d** macArray;
  Box2d arrayDelimiter;
  Box2d inflow(0,0,1,ny-2);
  PetscScalar ux,uy;

  ierr = solver.getInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);

  for (auto jj : arrayDelimiter.yRange){
    for (auto ii : arrayDelimiter.xRange){
      vf(0,ii,jj,ux,uy);
      macArray[jj][ii].ux = ux;
      macArray[jj][ii].uy = uy;
      // macArray[jj][ii].rho = rho0 + deltaRho;
    }
  }

  ierr = solver.restoreInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);

  ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
  ierr = solver.collideAndSwap(); CHKERRQ(ierr);

  AdjointLBSolver adjSolver;
  // ierr = solver.makeCodiAdjoint(&adjSolver); CHKERRQ(ierr);
  ierr = solver.makeSourceAdjoint(&adjSolver); CHKERRQ(ierr);

  /* Objective */
  PressureDrop2d<CollisionOperator> pres(inflow,Box2d(nx-1,nx-1,th+1,2*th),
                                         solver);
  /* Finite difference check */
  PetscInt numSteps = 1000;
  ComputeObjective obj(numSteps);
  NaiveCheckpointing check(numSteps+1,solver.distributionsLocal);
  // LBStaticCheckpointing check(100,numSteps+1,solver);
  // DensityFilter2d filter(7,Box2d(30,40,1,ny-2),solver,adjSolver);
  VariableDoubleProjectionFilter2d filter(4,2,Box2d(30,40,1,ny-2),solver,adjSolver);
  filter.setFirstBeta(8.);
  filter.setSecondBeta(4.);
  filter.setFirstEta(0.5);
  ierr = filter.setSecondEtaValuesFromFunction(EtaComputation()); CHKERRQ(ierr);
  PetscScalar value;
  obj(solver,pres,&value,true);
  PetscPrintf(PETSC_COMM_WORLD,"Pressure drop: %f\n",value);
  // LBDynamicCheckpointing check(500,solver);
  ComputeObjectiveAndSens objAndSens(check,numSteps);
  FiniteDifferenceCheck2d fd(obj,objAndSens,1e-5);
  // fd.checkDomain(solver,adjSolver,pres,Box2d(nx-10,nx-1,th+1,2*th),&filter);
  fd.checkDomain(solver,adjSolver,pres,Box2d(50,100,33,66));

  /* Main loop */

  // double t1,t2;
  // PetscInt step = 0;
  // solver.outputMacros();
  // PetscScalar totalDrop = 0.;
  // PetscScalar tstepDrop;
  // PetscInt numSteps = 1000;
  // PetscScalar timeAveraging = (PetscScalar) 1./numSteps;
  // NaiveCheckpointing check(numSteps+1,solver.distributionsLocal);
  // // LBStaticCheckpointing check(100,numSteps+1,solver);
  // check.makeCheckpoint(step,solver.distributionsLocal,nullptr);

  // t1 = MPI_Wtime();

  // while (step < numSteps){
  //   ierr = solver.streamBySwapping(); CHKERRQ(ierr);
  //   ++step;
  //   ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
  //   ierr = solver.collideAndSwap(); CHKERRQ(ierr);
  //   ierr = solver.computeMacroFunctional(pres,timeAveraging,&tstepDrop); CHKERRQ(ierr);
  //   totalDrop += tstepDrop;

  //   if (step % 50 == 0){
  //     ierr = solver.outputMacros(); CHKERRQ(ierr);
  //   }
  // }

  // t2 = MPI_Wtime();
  // PetscPrintf(PETSC_COMM_WORLD,"Pressure drop: %f\n",totalDrop);
  // PetscPrintf(PETSC_COMM_WORLD,"Forward time: %f\n",t2-t1);

  // // Adjoint
  // t1 = MPI_Wtime();
  // while (step > 0){
  //   ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
  //   ierr = adjSolver.adjointCollide(solver.distributionsLocal,
  //                                   solver.materialLocal,
  //                                   pres,1./numSteps); CHKERRQ(ierr);
  //   ierr = adjSolver.adjointStream(pres);
  //   --step;
  // }
  // t2 = MPI_Wtime();
  // PetscPrintf(PETSC_COMM_WORLD,"Adjoint time: %f\n",t2-t1);
}
