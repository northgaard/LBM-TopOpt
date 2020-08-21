#include "TOpTEnLB2d.hh"
#include "external/SobolGenerator.hh"

static char help[] = "The regenerator heat problem.\n";
static char outputFolder[] = "/home/sebnorg/WorkData/RobustRegenerator2d/";

using IsoLattice = D2Q9;
using ThermalLattice = D2Q4;
// using BaseIsoCollision = IncompressibleBGK2d<IsoLattice>;
// using BaseIsoCollision = StandardBGK2d<IsoLattice>;
using BaseIsoCollision = IncompressibleMRT2d<IsoLattice>;
// using BaseIsoCollision = Cascaded2d<IsoLattice>;
using IsoInterpolation = BPInterpolation<>;
using IsoCollision = PartialBouncebackCollision<BaseIsoCollision,IsoInterpolation>;
using BaseThermalCollision = ThermalBGK2d<ThermalLattice>;
using ThermalInterpolation = RampInterpolation<>;
using ThermalCollision = ThermalObstacleCollision<BaseThermalCollision,ThermalInterpolation>;
using CollisionOperator = ThermalDualCollision<IsoCollision,ThermalCollision>;

// Dirty dirty globals
static PetscInt stepsToSteady;

struct AEtaComputation {

  PetscScalar operator()(PetscInt idx, PetscInt idy) const
  {
    PetscScalar wave = 0.;
    auto inv = static_cast<PetscScalar>(numXWaves+numYWaves);
    inv = 1./sqrt(inv);
    for (PetscInt ii = 0; ii < numXWaves; ++ii){
      wave += inv*(coef[2*ii]*std::sin(2.*M_PI*omega_x[ii]*idx/nx) +
                   coef[2*ii+1]*std::cos(2.*M_PI*omega_x[ii]*idx/nx));
    }
    for (PetscInt ii = numXWaves; ii < numXWaves+numYWaves; ++ii){
      wave += inv*(coef[2*ii]*std::sin(2.*M_PI*omega_y[ii-numXWaves]*idy/ny) +
                   coef[2*ii+1]*std::cos(2.*M_PI*omega_y[ii-numXWaves]*idy/ny));
    }
    return (etaMax - etaMin)*0.5*(1. + std::erf(wave/std::sqrt(2.))) + etaMin;
  }
  PetscReal* coef;
  PetscReal* omega_x;
  PetscReal* omega_y;
  PetscInt nx;
  PetscInt ny;
  PetscReal etaMin, etaMax;
  PetscInt numXWaves;
  PetscInt numYWaves;
};

struct EtaComputation {

  PetscScalar operator()(PetscInt idx, PetscInt idy) const
  {
    PetscScalar wave = 1./std::sqrt(2.)*(coef[0]*std::sin(2.*M_PI*omega_y*idy/ny)
                                         + coef[1]*std::cos(2.*M_PI*omega_y*idy/ny))
      + 1./std::sqrt(2.)*(coef[2]*std::sin(2.*M_PI*omega_x*idx/nx)
                          + coef[3]*std::cos(2.*M_PI*omega_x*idx/nx));
    return (etaMax - etaMin)*0.5*(1. + std::erf(wave/std::sqrt(2.))) + etaMin;
  }
  PetscReal* coef;
  PetscReal omega_y;
  PetscReal omega_x;
  PetscInt nx;
  PetscInt ny;
  PetscReal etaMin, etaMax;
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

  void operator()(PetscInt,PetscInt,PetscInt,PetscScalar& rho) const
  {
    rho = rho0;
  }
  PetscScalar rho0;
};

struct StageTemp {

  void operator()(PetscInt t, PetscInt, PetscInt, PetscScalar& T) const
  {
    // Still computing steady state
    if (t <= stepsToSteady){
      T = T0;
    } else { // Ramp up temperature
      PetscInt t_u = t-stepsToSteady;
      PetscScalar mult;
      if (t_u < rampUpTime){
        mult = sin(M_PI*t_u/(2.*rampUpTime));
      } else {
        mult = 1.;
      }
      T = T0 + mult*deltaT;
    }
  }
  PetscScalar T0,deltaT;
  PetscInt rampUpTime;
};

class ComputeSteadyState {

public:
  ComputeSteadyState(Vec _v)
  {
    VecDuplicate(_v,&prev);
    VecDuplicate(_v,&diffVec);
  }
  PetscErrorCode operator()(NewLBSolver2d& solver)
  {
    PetscErrorCode ierr;
    PetscScalar tol = 1e-7;
    PetscScalar norm = 1.;
    PetscScalar procNorm;
    PetscInt step = 0;
    static constexpr PetscInt chkFreq = 100;
    PetscFunctionBeginUser;
    ierr = solver.initializeAtEquilibrium(); CHKERRQ(ierr);
    ierr = solver.setCurrentTimestep(0); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);

    // Compute steady state
    while (true){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      if (step % chkFreq == 0){
        ierr = VecWAXPY(diffVec,-1.,solver.distributionsLocal,prev); CHKERRQ(ierr);
        ierr = VecNorm(diffVec,NORM_INFINITY,&procNorm); CHKERRQ(ierr);
        MPI_Allreduce(&procNorm,&norm,1,MPIU_SCALAR,MPI_MAX,PETSC_COMM_WORLD);
        if (norm < tol){
          break;
        }
      }
      if (step % chkFreq == chkFreq-1){
        ierr = VecCopy(solver.distributionsLocal,prev); CHKERRQ(ierr);
      }
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    }
    PetscPrintf(PETSC_COMM_WORLD,"Number of steps to steady state: %d\n",step);
    PetscFunctionReturn(0);
  }
private:
  Vec prev,diffVec;
};

class ComputeObjective {

public:
  ComputeObjective(PetscInt _nt, Vec _steady) : numUnsteadySteps(_nt)
  {
    timeAveraging = (PetscScalar) 1./numUnsteadySteps;
    VecDuplicate(_steady,&lastSteady);
    VecCopy(_steady,lastSteady);
    VecDuplicate(lastSteady,&prev);
    VecDuplicate(lastSteady,&diffVec);
  }
  PetscErrorCode setSteadyInit(Vec init)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = VecCopy(init,lastSteady); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()(NewLBSolver2d& solver, const LBMacroFunctional& func,
                            const LBMacroFunctional& effFunc,
                            PetscScalar* value, PetscScalar* effValue,
                            bool outputMacros = false)
  {
    PetscErrorCode ierr;
    PetscScalar tol = 1e-6;
    PetscScalar norm = 1.;
    PetscScalar procNorm;
    PetscInt step = 0;
    static constexpr PetscInt chkFreq = 100;
    PetscFunctionBeginUser;
    stepsToSteady = std::numeric_limits<PetscInt>::max();
    ierr = solver.setDistributions(lastSteady); CHKERRQ(ierr);
    ierr = solver.setCurrentTimestep(0); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);

    // Compute steady state
    while (true){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      if (step % chkFreq == 0){
        ierr = VecWAXPY(diffVec,-1.,solver.distributionsLocal,prev); CHKERRQ(ierr);
        ierr = VecNorm(diffVec,NORM_INFINITY,&procNorm); CHKERRQ(ierr);
        MPI_Allreduce(&procNorm,&norm,1,MPIU_SCALAR,MPI_MAX,PETSC_COMM_WORLD);
        PetscPrintf(PETSC_COMM_WORLD,"Norm difference at step %d: %f\n",step,norm);
        if (norm < tol){
          ierr = VecCopy(solver.distributionsLocal,lastSteady); CHKERRQ(ierr);
          if (outputMacros){
            char steadyOut[50];
            sprintf(steadyOut,"Steady_data_step_%d.vts",step);
            ierr = solver.outputMacros(outputFolder,steadyOut); CHKERRQ(ierr);
          }
          break;
        }
      }
      if (step % chkFreq == chkFreq-1){
        ierr = VecCopy(solver.distributionsLocal,prev); CHKERRQ(ierr);
      }
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    }
    // Compute unsteady
    stepsToSteady = step;
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    PetscInt endStep = step + numUnsteadySteps;
    PetscScalar tvalue;
    *value = 0.;
    *effValue = 0.;
    while (step < endStep){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
      if (outputMacros && (step-stepsToSteady) % 200 == 0){
        char output[50];
        sprintf(output,"Macros_unsteady_step_%d.vts",step-stepsToSteady);
        ierr = solver.outputMacros(outputFolder,output); CHKERRQ(ierr);
      }
      ierr = solver.computeMacroFunctional(func,timeAveraging,&tvalue);
      *value += tvalue;
    }
    ierr = solver.computeMacroFunctional(effFunc,1.,effValue); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  Vec lastSteady;
private:
  Vec prev,diffVec;
  PetscInt numUnsteadySteps;
  PetscScalar timeAveraging;
};

class ComputeObjectiveAndSens {

public:
  ComputeObjectiveAndSens(LBDynamicCheckpointing& _chk, PetscInt _nt, Vec _v)
    : check(_chk), numUnsteadySteps(_nt)
  {
    timeAveraging = (PetscScalar) 1./numUnsteadySteps;
    VecDuplicate(_v,&lastSteady);
    VecCopy(_v,lastSteady);
    VecDuplicate(lastSteady,&prev);
    VecDuplicate(lastSteady,&diffVec);
  }
  PetscErrorCode setSteadyInit(Vec init)
  {
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    ierr = VecCopy(init,lastSteady); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()(NewLBSolver2d& solver, AdjointLBSolver& adjSolver,
                            AdjointLBSolver& conAdjSolver, const LBMacroFunctional& func,
                            const LBMacroFunctional& conFunc, PetscScalar* value,
                            PetscScalar* conValue, Vec* sens, Vec* conSens)
  {
    PetscErrorCode ierr;
    PetscScalar tol = 1e-7;
    PetscScalar norm;
    PetscScalar procNorm;
    PetscInt step = 0;
    double t1,t2;
    static constexpr PetscInt chkFreq = 100;
    static constexpr PetscInt maxSteps = 150000;
    PetscFunctionBeginUser;
    stepsToSteady = std::numeric_limits<PetscInt>::max();
    check.reset(nullptr);
    ierr = solver.setDistributions(lastSteady); CHKERRQ(ierr);
    ierr = solver.setCurrentTimestep(0); CHKERRQ(ierr);
    ierr = check.makeCheckpoint(0,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);

    /*
      Primal solve
    */

    t1 = MPI_Wtime();

    // Compute steady state
    while (true){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      if (step % chkFreq == 0){
        ierr = VecWAXPY(diffVec,-1.,solver.distributionsLocal,prev); CHKERRQ(ierr);
        ierr = VecNorm(diffVec,NORM_INFINITY,&procNorm); CHKERRQ(ierr);
        MPI_Allreduce(&procNorm,&norm,1,MPIU_SCALAR,MPI_MAX,PETSC_COMM_WORLD);
        if (norm < tol || step == maxSteps){
          ierr = VecCopy(solver.distributionsLocal,lastSteady); CHKERRQ(ierr);
          break;
        }
      }
      if (step % chkFreq == chkFreq-1){
        ierr = VecCopy(solver.distributionsLocal,prev); CHKERRQ(ierr);
      }
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    }

    // Compute unsteady
    stepsToSteady = step;
    if (stepsToSteady == maxSteps){
      PetscPrintf(PETSC_COMM_WORLD,"Reached max number of steps for steady computation, norm difference is %f\n",norm);
    }
    PetscPrintf(PETSC_COMM_WORLD,"Number of steps to steady state: %d\n",stepsToSteady);
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    PetscInt endStep = stepsToSteady + numUnsteadySteps;
    PetscScalar tvalue, ctvalue;
    *value = 0.;
    *conValue = 0.;
    while (step < endStep){
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = check.makeCheckpoint(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
      // ierr = solver.computeMacroFunctional(func,timeAveraging,&tvalue); CHKERRQ(ierr);
      ierr = solver.computeMacroFunctional(conFunc,timeAveraging,&ctvalue); CHKERRQ(ierr);
      // *value += tvalue;
      *conValue += ctvalue;
    }
    ierr = solver.computeMacroFunctional(func,1.,&tvalue); CHKERRQ(ierr);
    *value = 1./tvalue;

    t2 = MPI_Wtime();

    PetscPrintf(PETSC_COMM_WORLD,"Forward time: %f\n",t2-t1);

    /*
      Adjoint solve
    */
    adjSolver.resetAdjoints();
    adjSolver.setCurrentTimestep(step);
    conAdjSolver.resetAdjoints();
    conAdjSolver.setCurrentTimestep(step);

    t1 = MPI_Wtime();

    // Unsteady part
    while (step > stepsToSteady){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal);
      CHKERRQ(ierr);
      if (step == endStep){
        ierr = adjSolver.computeCollideSource(func,1.,solver.distributionsLocal,
                                              solver.materialLocal); CHKERRQ(ierr);
        ierr = adjSolver.computeSensitivitySource(func,1.,solver.distributionsLocal,
                                                  solver.materialLocal); CHKERRQ(ierr);
      }
      ierr = adjSolver.adjointStream(solver.distributionsLocal); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal);
      CHKERRQ(ierr);
      ierr = conAdjSolver.computeCollideSource(conFunc,timeAveraging,solver.distributionsLocal,
                                               solver.materialLocal); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointStream(solver.distributionsLocal); CHKERRQ(ierr);
      --step;
    }

    // Steady part
    while (step > 0){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal);
      CHKERRQ(ierr);
      ierr = adjSolver.adjointStream(solver.distributionsLocal); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal);
      CHKERRQ(ierr);
      ierr = conAdjSolver.adjointStream(solver.distributionsLocal); CHKERRQ(ierr);
      --step;
    }

    t2 = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Adjoint time: %f\n",t2-t1);
    adjSolver.getSensitivities(sens);
    // Rescale objective sensitivities since we are doing 1/efficiency
    ierr = VecScale(*sens,-1./(tvalue*tvalue)); CHKERRQ(ierr);
    conAdjSolver.getSensitivities(conSens);
    PetscFunctionReturn(0);
  }
  Vec lastSteady;
private:
  LBDynamicCheckpointing& check;
  PetscInt numUnsteadySteps;
  PetscScalar timeAveraging;
  Vec prev,diffVec;
};

int main(int argc, char *argv[])
{
  TopTenInit init(argc,argv,help);
  PetscErrorCode ierr;
  constexpr PetscMPIInt numCases = 1;
  ThermalLBSolver2d<CollisionOperator> solver;

  MPIEvenSplit split(numCases);
  MPI_Comm splitComm;
  PetscMPIInt splitID;
  PetscMPIInt* roots;
  split.getCommunicationInfo(&splitComm,&splitID,&roots);

  NewLBSolverInfo2d info;
  // PetscInt nx = 550;
  PetscInt nx = 350;
  PetscInt ny = 140;
  info.nx = nx;
  info.ny = ny;

  IncompressibleFlowParameters isoPar;
  isoPar.velocityChar = 0.01;
  isoPar.ReynoldsNumber = 30;
  isoPar.lengthChar = ny;
  ThermalFlowParameters thermalPar;
  thermalPar.incPar = isoPar;
  thermalPar.PrandtlNumber = 0.71;

  BaseIsoCollision baseIsoOp(isoPar);
  IsoInterpolation isoInterp(1.);
  IsoCollision isoOp(baseIsoOp,isoInterp);

  BaseThermalCollision baseThermalOp(thermalPar);
  // Water/Aluminium Ck approx 0.0025
  // PetscScalar Ck = 0.01;
  PetscScalar Ck = 0.0025;
  ThermalInterpolation thermalInterp(1.,Ck);
  ThermalCollision thermalOp(baseThermalOp,thermalInterp);
  CollisionOperator colOp(isoOp,thermalOp);

  solver.make(splitComm,info,colOp);

  // Isothermal boundaries
  OnGridBounceback2d<IsoLattice> bb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),bb); CHKERRQ(ierr);

  // InflowVelocity vf;
  // vf.vmax = isoPar.velocityChar;
  // vf.mid = (PetscScalar) (ny-1)/2.;
  // ZouHeVelocity2d<IsoCollision,IsoLattice,InflowVelocity> zhv(vf);
  // ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),zhv); CHKERRQ(ierr);

  PetscScalar rho0 = 1.;
  ConstPres outP;
  outP.rho0 = rho0;
  ZouHePressure2d<IsoLattice,IsoCollision::Equilibrium,ConstPres> zhp(outP);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,1,ny-2),zhp); CHKERRQ(ierr);
  // Prescibed pressure drop
  PetscScalar deltaRho = 0.05;
  PetscScalar rhoIn = rho0 + deltaRho;
  ConstPres inP;
  inP.rho0 = rhoIn;
  ZouHePressure2d<IsoLattice,IsoCollision::Equilibrium,ConstPres> zhpIn(inP);
  ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),zhpIn); CHKERRQ(ierr);

  // Thermal boundaries
  ThermalBounceback2d<IsoLattice,ThermalLattice> tbb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),tbb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),tbb); CHKERRQ(ierr);

  PetscScalar baseTemp = 4.;
  PetscScalar dT = -2.;
  StageTemp st;
  st.T0 = baseTemp;
  st.deltaT = dT;
  st.rampUpTime = 200;
  ThermalZouHe2d<IsoLattice,ThermalLattice,StageTemp> stzh(st);
  ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),stzh); CHKERRQ(ierr);

  ThermalZouHeNeumann2d<IsoLattice,ThermalLattice> tzhn;
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,1,ny-2),tzhn); CHKERRQ(ierr);

  ThermalMacros2d mac;
  mac.rho = 1.;
  mac.ux = 0.;
  mac.uy = 0.;
  mac.T = baseTemp;
  ierr = solver.uniformMacroInitialization(mac); CHKERRQ(ierr);

  ThermalMacros2d** macArray;
  Box2d arrayDelimiter;
  Box2d inflow(0,0,1,ny-2);
  // PetscScalar ux,uy;

  // Velocity inflow
  // ierr = solver.getInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);
  // for (auto jj : arrayDelimiter.yRange){
  //   for (auto ii : arrayDelimiter.xRange){
  //     vf(0,ii,jj,ux,uy);
  //     macArray[jj][ii].ux = ux;
  //     macArray[jj][ii].uy = uy;
  //   }
  // }
  // ierr = solver.restoreInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);

  // Pressure
  ierr = solver.getInitialMacroArray(solver.getBoundingBox(),&arrayDelimiter,
                                     &macArray); CHKERRQ(ierr);
  for (auto jj : arrayDelimiter.yRange){
    for (auto ii : arrayDelimiter.xRange){
      macArray[jj][ii].rho = (rho0 - rhoIn)*ii/(nx-1) + rhoIn;
    }
  }
  ierr = solver.restoreInitialMacroArray(solver.getBoundingBox(),&arrayDelimiter,
                                         &macArray); CHKERRQ(ierr);

  PetscScalar solid[1] = {0.};
  PetscInt freeLen = 50;
  PetscInt oFreeLen = 100;
  // PetscInt plateSep = ny/7;
  // Initial plates (ny = 140)
  // ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,0,plateSep-1),solid);
  // CHKERRQ(ierr);
  // ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,2*plateSep,3*plateSep-1),solid);
  // CHKERRQ(ierr);
  // ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,4*plateSep,5*plateSep-1),solid);
  // CHKERRQ(ierr);
  // ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,6*plateSep,ny-1),solid); CHKERRQ(ierr);
  PetscInt plateSep = ny/14;
  PetscInt hps = plateSep/2;
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,0,hps-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+plateSep,hps+2*plateSep-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+3*plateSep,hps+4*plateSep-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+5*plateSep,hps+6*plateSep-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+7*plateSep,hps+8*plateSep-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+9*plateSep,hps+10*plateSep-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+11*plateSep,hps+12*plateSep-1),solid);
  CHKERRQ(ierr);
  ierr = solver.setUniformFieldValues(Box2d(freeLen,nx-1-oFreeLen,hps+13*plateSep,ny-1),solid);
  CHKERRQ(ierr);

  ierr = solver.outputFields(outputFolder,"Base.vts"); CHKERRQ(ierr);

  // Set objective
  // PetscInt xStart = nx-oFreeLen + 40;
  // PetscInt domLength = 10;
  // DomainTemperature2d<CollisionOperator> objective(Box2d(xStart,xStart+domLength-1,0,ny-1),
  //                                                  solver.getLocalBoundingBox());
  RegeneratorEfficiency2d<CollisionOperator> efficiency(Box2d(freeLen,nx-1-oFreeLen,0,ny-1),
                                                        solver,baseTemp,baseTemp+dT);
  ThermalOutflowEast2d<CollisionOperator> tflow(Box2d(nx-1,nx-1,1,ny-2),solver);
  // Pressure drop constraint
  // PressureDrop2d<IsoCollision> stateConstraint(Box2d(0,0,1,ny-2),Box2d(nx-1,nx-1,1,ny-2),
  //                                              solver);

  PetscInt unsteadySteps = 35000;
  // LBDynamicCheckpointing check(2000,solver);
  solver.initializeAtEquilibrium();
  // ComputeObjectiveAndSens opt(check,unsteadySteps,solver.distributionsLocal);
  ComputeObjective objComp(unsteadySteps,solver.distributionsLocal);
  PetscScalar tempValue;
  PetscScalar effValue;
  // // Compute reference pressure
  objComp(solver,tflow,efficiency,&tempValue,&effValue,false);
  PetscPrintf(PETSC_COMM_WORLD,"Thermal outflow: %f\nEfficiency: %f\n",tempValue,effValue);
  // // Set steady state for optimization (since we just computed it anyway)
  // ierr = opt.setSteadyInit(objComp.lastSteady); CHKERRQ(ierr);
  // PetscPrintf(PETSC_COMM_WORLD,"Efficiency for initial plates: %f\n",tempValue);

  AdjointLBSolver adjSolver;
  ierr = solver.makeSourceAdjoint(&adjSolver); CHKERRQ(ierr);
  AdjointLBSolver presAdjSolver;
  ierr = solver.makeSourceAdjoint(&presAdjSolver); CHKERRQ(ierr);
  Box2d designDomainBox(freeLen,nx-1-oFreeLen,0,ny-1);

  AEtaComputation etaComp;
  SobolSeedType seed = 10000 + splitID;
  SobolGenerator gen(seed);
  constexpr PetscInt numWaves_x = 12;
  constexpr PetscInt numWaves_y = 12;
  auto coef = gen.getNormalDistributedSequence(2*numWaves_x+2*numWaves_y,0.,1.);
  etaComp.coef = &coef[0];
  // PetscReal omega_x[numWaves] = {2.,4.,8.,16.,32.,64.};
  // PetscReal omega_y[numWaves] = {2.,4.,8.,16.,32.,64.};
  // PetscReal omega_x[numWaves] = {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.};
  // PetscReal omega_y[numWaves] = {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.};
  PetscReal omega_x[numWaves_x] = {2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.};
  PetscReal omega_y[numWaves_y] = {2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.};
  // PetscReal omega_x[numWaves] = {2.,4.,6.,8.,10.,20.};
  // PetscReal omega_y[numWaves] = {2.,4.,6.,8.,10.,20.};
  // PetscReal omega_x[3] = {5.,10.,15.};
  // PetscReal omega_x[numWaves_x] = {1.,2.,3.};
  // PetscReal omega_y[numWaves_y] = {0.1,0.2,0.3,0.4,0.5,1.,1.5,2.,2.5,3.,5.};
  etaComp.omega_x = &omega_x[0];
  etaComp.omega_y = &omega_y[0];
  etaComp.nx = designDomainBox.getNx();
  etaComp.ny = designDomainBox.getNy();
  etaComp.numXWaves = numWaves_x;
  etaComp.numYWaves = numWaves_y;
  // etaComp.etaMin = 0.05;
  // etaComp.etaMax = 0.95;
  etaComp.etaMin = 0.2;
  etaComp.etaMax = 0.8;
  // EtaComputation etaComp;
  // SobolSeedType seed = 10000 + splitID;
  // SobolGenerator gen(seed);
  // auto coef = gen.getNormalDistributedSequence(4,0.,1.);
  // etaComp.coef = &coef[0];
  // etaComp.omega_y = 64;
  // etaComp.omega_x = 16;
  // etaComp.nx = designDomainBox.getNx();
  // etaComp.ny = designDomainBox.getNy();
  // etaComp.etaMin = 0.3;
  // etaComp.etaMax = 0.7;
  // PetscPrintf(splitComm,"Coefficients: %f %f %f %f\n",coef[0],coef[1],coef[2],coef[3]);

  Vec fullDomain, designDomain;
  Vec physicalFullDomain, physicalDesignDomain;
  Vec sens, designSens, presSens, designPresSens;
  Vec vcSens, designVcSens;
  FullToDesignMapping2d map(designDomainBox,solver);
  // ProjectionFilter2d filter(4,Box2d(freeLen,nx-1-oFreeLen,0,ny-1),solver,adjSolver);
  // filter.setBeta(1.);
  // filter.setEta(0.5);
  VariableProjectionFilter2d filter(4,designDomainBox,solver,adjSolver);
  PetscScalar betaMax = std::pow(1.5,12.);
  PetscPrintf(PETSC_COMM_WORLD,"Beta max: %f\n",betaMax);
  filter.setBeta(betaMax);
  ierr = filter.setEtaValuesFromFunction(etaComp); CHKERRQ(ierr);
  // PetscScalar meanEta;
  // ierr = filter.meanEtaValue(&meanEta); CHKERRQ(ierr);
  // PetscPrintf(splitComm,"Mean eta, subcomm %d: %f\n",splitID,meanEta);

  // char output[50];
  // sprintf(output,"Eta_values_realization_%d.vts",splitID);
  // ierr = filter.outputEtaValues(outputFolder,output); CHKERRQ(ierr);

  ierr = solver.createDesignVec(&fullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&physicalFullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&vcSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designDomain); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designVcSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designPresSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&physicalDesignDomain); CHKERRQ(ierr);

  // ierr = map.mapFullToDesign(fullDomain,designDomain); CHKERRQ(ierr);

  // ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
  // ierr = filter.filterDesign(fullDomain,physicalFullDomain); CHKERRQ(ierr);
  // ierr = solver.setFieldFromDesignVec(physicalFullDomain); CHKERRQ(ierr);

  // sprintf(output,"Filtered_design_realization_%d.vts",splitID);
  // ierr = solver.outputFields(outputFolder,output);

  // VolumeConstraint vc(0.5,ConstrainFluid,designDomain); CHKERRQ(ierr);
  // PetscScalar volumeValue;
  // PetscScalar objValue;
  // PetscScalar presValue;
  // PetscScalar presConst = 1.2;
  // PetscInt numConstraints = 2;
  // PetscInt numDesign;
  // ierr = VecGetSize(designDomain,&numDesign); CHKERRQ(ierr);
  // MMA mma(numDesign,numConstraints,designDomain);
  // PetscScalar tol = 1e-4;
  // PetscScalar ch = 1.;
  // PetscInt maxIter = 1200;
  // PetscInt curIter = 0;
  // PetscScalar movlim = 0.1;
  // PetscInt incFreq = 50;
  // PetscScalar betaInc = 1.5;
  // char output[50];
  // Vec conSens[2];
  // Vec conDesignSens[2];
  // adjSolver.getSensitivities(&sens);
  // presAdjSolver.getSensitivities(&presSens);
  // conSens[0] = presSens;
  // conSens[1] = vcSens;
  // conDesignSens[0] = designPresSens;
  // conDesignSens[1] = designVcSens;
  // PetscScalar conValue[2];
  // PetscPrintf(PETSC_COMM_WORLD,"-----------------------\n");

  // // Optimization loop
  // while (ch > tol && curIter < maxIter){

  //   ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
  //   ierr = filter.filterDesign(fullDomain,physicalFullDomain); CHKERRQ(ierr);
  //   ierr = solver.setFieldFromDesignVec(physicalFullDomain); CHKERRQ(ierr);

  //   if (curIter % 50 == 0){
  //     sprintf(output,"Design_iter_%d.vts",curIter);
  //     ierr = solver.outputFields(outputFolder,output); CHKERRQ(ierr);
  //   }

  //   // Compute objective and sensitivities
  //   ierr = opt(solver,adjSolver,presAdjSolver,objective,stateConstraint,&objValue,&presValue,
  //              &sens,&presSens); CHKERRQ(ierr);
  //   // Rescale pressure drop constraint
  //   presValue = presValue/(presConst*refPressure) - 1.;
  //   ierr = VecScale(presSens,1./(presConst*refPressure)); CHKERRQ(ierr);

  //   // Volume constraint
  //   ierr = map.mapFullToDesign(physicalFullDomain,physicalDesignDomain); CHKERRQ(ierr);
  //   vc.computeFunctional(physicalDesignDomain,&volumeValue);
  //   vc.computeSensitivities(vcSens);

  //   // Filter sensitivities
  //   ierr = filter.filterSensitivities(sens,conSens,numConstraints); CHKERRQ(ierr);
  //   // Map to design domain
  //   ierr = map.mapFullToDesign(sens,designSens); CHKERRQ(ierr);
  //   ierr = map.mapFullToDesign(presSens,designPresSens); CHKERRQ(ierr);
  //   ierr = map.mapFullToDesign(vcSens,designVcSens); CHKERRQ(ierr);

  //   PetscPrintf(PETSC_COMM_WORLD,"Iter: %d, objective: %f, pressure constraint: %f, volume constraint: %f\n",
  //               curIter,objValue,presValue,volumeValue);
  //   conValue[0] = presValue;
  //   conValue[1] = volumeValue;
  //   mma.Update(designDomain,designSens,conValue,conDesignSens,movlim);
  //   mma.DesignChange(designDomain,ch);
  //   ++curIter;
  //   if (curIter % incFreq == 0 && curIter <= 600){
  //     filter.increaseBeta(betaInc);
  //   }
  //   PetscPrintf(PETSC_COMM_WORLD,"Design change: %f\n",ch);
  //   PetscPrintf(PETSC_COMM_WORLD,"-----------------------\n");
  // }

  // // Output final design
  // ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
  // ierr = filter.filterDesign(fullDomain,physicalFullDomain); CHKERRQ(ierr);
  // ierr = solver.setFieldFromDesignVec(physicalFullDomain); CHKERRQ(ierr);
  // ierr = solver.outputFields(outputFolder,"Final_design.vts"); CHKERRQ(ierr);

  // // Compute final objective/constraint and output macros
  // ierr = map.mapFullToDesign(physicalFullDomain,physicalDesignDomain); CHKERRQ(ierr);
  // vc.computeFunctional(physicalDesignDomain,&volumeValue);
  // ierr = objComp.setSteadyInit(opt.lastSteady); CHKERRQ(ierr);
  // objComp(solver,objective,stateConstraint,&objValue,&presValue,true);
  // presValue = presValue/(presConst*refPressure) - 1.;
  // PetscPrintf(PETSC_COMM_WORLD,"Final design: objective: %f, pressure constraint: %f, volume constraint: %f\n",objValue,presValue,volumeValue);

  ierr = VecDestroy(&fullDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&physicalFullDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&vcSens); CHKERRQ(ierr);
  ierr = VecDestroy(&designDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&designSens); CHKERRQ(ierr);
  ierr = VecDestroy(&designVcSens); CHKERRQ(ierr);
  ierr = VecDestroy(&physicalDesignDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&designPresSens); CHKERRQ(ierr);
}
