#include "TOpTEnLB2d.hh"

static char help[] = "The regenerator heat problem.\n";
static char outputFolder[] = "/home/sebnorg/WorkData/Regenerator/";

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

// Hack the adjoint implementation, this also needs to change!
class EmptyFunctional : public LBMacroFunctional {

public:
  void evaluate(PetscInt,const PetscScalar,void*,PetscScalar*) const override {}
  void adjointCollideSource(PetscInt,const PetscScalar,void*,void*) const override {}
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
    PetscScalar tol = 1e-6;
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
                            const LBMacroFunctional& conFunc,
                            PetscScalar* value, PetscScalar* conValue,
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
        // PetscPrintf(PETSC_COMM_WORLD,"Norm difference at step %d: %f\n",step,norm);
        // char steadyOut[50];
        // sprintf(steadyOut,"Steady_data_step_%d.vts",step);
        // ierr = solver.outputMacros(outputFolder,steadyOut); CHKERRQ(ierr);
        if (norm < tol){
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
    ierr = solver.collideAndSwap(); CHKERRQ(ierr);
    PetscInt endStep = step + numUnsteadySteps;
    PetscScalar tvalue;
    PetscScalar ctvalue;
    *value = 0.;
    *conValue = 0.;
    while (step < endStep){
      if (outputMacros && (step-stepsToSteady) % 200 == 0){
        char output[50];
        sprintf(output,"Macros_unsteady_step_%d.vts",step-stepsToSteady);
        ierr = solver.outputMacros(outputFolder,output); CHKERRQ(ierr);
      }
      ierr = solver.streamBySwapping(); CHKERRQ(ierr);
      ++step;
      ierr = solver.collideAndSwap(); CHKERRQ(ierr);
      ierr = solver.computeMacroFunctional(func,timeAveraging,&tvalue);
      ierr = solver.computeMacroFunctional(conFunc,timeAveraging,&ctvalue);
      *value += tvalue;
      *conValue += ctvalue;
    }
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
    PetscScalar tol = 1e-6;
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
      ierr = solver.computeMacroFunctional(func,timeAveraging,&tvalue); CHKERRQ(ierr);
      ierr = solver.computeMacroFunctional(conFunc,timeAveraging,&ctvalue); CHKERRQ(ierr);
      *value += tvalue;
      *conValue += ctvalue;
    }

    t2 = MPI_Wtime();

    PetscPrintf(PETSC_COMM_WORLD,"Forward time: %f\n",t2-t1);

    /*
      Adjoint solve
    */
    adjSolver.resetAdjoints();
    adjSolver.setCurrentTimestep(step);
    conAdjSolver.resetAdjoints();
    conAdjSolver.setCurrentTimestep(step);
    EmptyFunctional empt;

    t1 = MPI_Wtime();

    // Unsteady part
    while (step > stepsToSteady){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal,
                                      func,timeAveraging); CHKERRQ(ierr);
      ierr = adjSolver.adjointStream(func); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal,
                                         conFunc,timeAveraging); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointStream(conFunc); CHKERRQ(ierr);
      --step;
    }

    // Steady part
    while (step > 0){
      ierr = check.getForwardSolution(step,solver.distributionsLocal,nullptr); CHKERRQ(ierr);
      ierr = adjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal,
                                      empt,timeAveraging); CHKERRQ(ierr);
      ierr = adjSolver.adjointStream(empt); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointCollide(solver.distributionsLocal,solver.materialLocal,
                                         empt,timeAveraging); CHKERRQ(ierr);
      ierr = conAdjSolver.adjointStream(empt); CHKERRQ(ierr);
      --step;
    }

    t2 = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Adjoint time: %f\n",t2-t1);
    adjSolver.getSensitivities(sens);
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
  ThermalLBSolver2d<CollisionOperator> solver;

  NewLBSolverInfo2d info;
  PetscInt nx = 550;
  PetscInt ny = 140;
  info.nx = nx;
  info.ny = ny;

  IncompressibleFlowParameters isoPar;
  isoPar.velocityChar = 0.1;
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

  solver.make(info,colOp);

  // Isothermal boundaries
  OnGridBounceback2d<IsoLattice> bb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),bb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),bb); CHKERRQ(ierr);

  InflowVelocity vf;
  vf.vmax = isoPar.velocityChar;
  vf.mid = (PetscScalar) (ny-1)/2.;
  ZouHeVelocity2d<IsoCollision,IsoLattice,InflowVelocity> zhv(vf);
  ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),zhv); CHKERRQ(ierr);

  PetscScalar rho0 = 1.;
  ConstPres outP;
  outP.rho0 = rho0;
  ZouHePressure2d<IsoLattice,IsoCollision::Equilibrium,ConstPres> zhp(outP);
  ierr = solver.addBoundaryCondition(Box2d(nx-1,nx-1,1,ny-2),zhp); CHKERRQ(ierr);
  // Prescibed pressure drop
  // PetscScalar deltaRho = 2.0;
  // PetscScalar rhoIn = rho0 + deltaRho;
  // ConstPres inP;
  // inP.rho0 = rhoIn;
  // ZouHePressure2d<IsoLattice,IsoCollision::Equilibrium,ConstPres> zhpIn(inP);
  // ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),zhpIn); CHKERRQ(ierr);

  // Thermal boundaries
  ThermalBounceback2d<IsoLattice,ThermalLattice> tbb;
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,0,0),tbb); CHKERRQ(ierr);
  ierr = solver.addBoundaryCondition(Box2d(0,nx-1,ny-1,ny-1),tbb); CHKERRQ(ierr);

  // ConstTemp ct;
  PetscScalar baseTemp = 0.;
  // ct.T0 = baseTemp;
  // ThermalZouHe2d<IsoLattice,ThermalLattice,ConstTemp> tzh(ct);
  // ierr = solver.addBoundaryCondition(Box2d(0,0,1,ny-2),tzh); CHKERRQ(ierr);
  // Unsteady solver cranks up the temperature!
  // RampUpTemp rtemp;
  // rtemp.T0 = 2.;
  // rtemp.rampUpTime = 200;
  // ThermalZouHe2d<IsoLattice,ThermalLattice,RampUpTemp> rtzh(rtemp);
  StageTemp st;
  st.T0 = baseTemp;
  st.deltaT = 2.;
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
  PetscScalar ux,uy;

  // Velocity inflow
  ierr = solver.getInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);
  for (auto jj : arrayDelimiter.yRange){
    for (auto ii : arrayDelimiter.xRange){
      vf(0,ii,jj,ux,uy);
      macArray[jj][ii].ux = ux;
      macArray[jj][ii].uy = uy;
    }
  }
  ierr = solver.restoreInitialMacroArray(inflow,&arrayDelimiter,&macArray); CHKERRQ(ierr);

  // Pressure
  // ierr = solver.getInitialMacroArray(solver.getBoundingBox(),&arrayDelimiter,
  //                                    &macArray); CHKERRQ(ierr);
  // for (auto jj : arrayDelimiter.yRange){
  //   for (auto ii : arrayDelimiter.xRange){
  //     macArray[jj][ii].rho = (rho0 - rhoIn)*ii/(nx-1) + rhoIn;
  //   }
  // }
  // ierr = solver.restoreInitialMacroArray(solver.getBoundingBox(),&arrayDelimiter,
  //                                        &macArray); CHKERRQ(ierr);

  // PetscScalar solid[1] = {0.};
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

  // Set objective
  PetscInt xStart = nx-oFreeLen + 40;
  PetscInt domLength = 10;
  DomainTemperature2d<CollisionOperator> objective(Box2d(xStart,xStart+domLength-1,0,ny-1),
                                                   solver.getLocalBoundingBox());
  // Pressure drop constraint
  PressureDrop2d<IsoCollision> stateConstraint(Box2d(0,0,1,ny-2),Box2d(nx-1,nx-1,1,ny-2),
                                               solver.getLocalBoundingBox());

  PetscInt unsteadySteps = 25000;
  LBDynamicCheckpointing check(50,solver);
  solver.initializeAtEquilibrium();
  ComputeObjectiveAndSens opt(check,unsteadySteps,solver.distributionsLocal);
  ComputeObjective objComp(unsteadySteps,solver.distributionsLocal);
  PetscScalar tempValue;
  PetscScalar refPressure;
  // Compute reference pressure
  objComp(solver,objective,stateConstraint,&tempValue,&refPressure);
  PetscPrintf(PETSC_COMM_WORLD,"Objective: %f, reference pressure: %f\n",tempValue,refPressure);
  // Set steady state for optimization (since we just computed it anyway)
  ierr = opt.setSteadyInit(objComp.lastSteady); CHKERRQ(ierr);

  AdjointLBSolver adjSolver;
  ierr = solver.makeSourceAdjoint(&adjSolver); CHKERRQ(ierr);
  AdjointLBSolver presAdjSolver;
  ierr = solver.makeSourceAdjoint(&presAdjSolver); CHKERRQ(ierr);

  Vec fullDomain, designDomain;
  Vec physicalFullDomain, physicalDesignDomain;
  Vec sens, designSens, presSens, designPresSens;
  Vec vcSens, designVcSens;
  FullToDesignMapping2d map(Box2d(freeLen,nx-1-oFreeLen,0,ny-1),solver);
  DensityFilter2d filter(3,Box2d(freeLen,nx-1-oFreeLen,0,ny-1),solver,adjSolver);

  ierr = solver.createDesignVec(&fullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&physicalFullDomain); CHKERRQ(ierr);
  ierr = solver.createDesignVec(&vcSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designDomain); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designVcSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&designPresSens); CHKERRQ(ierr);
  ierr = map.createDesignDomainVec(&physicalDesignDomain); CHKERRQ(ierr);

  ierr = map.mapFullToDesign(fullDomain,designDomain); CHKERRQ(ierr);

  VolumeConstraint vc(0.5,ConstrainFluid,designDomain); CHKERRQ(ierr);
  PetscScalar volumeValue;
  PetscScalar objValue;
  PetscScalar presValue;
  PetscScalar presConst = 0.8;
  PetscInt numConstraints = 2;
  PetscInt numDesign;
  ierr = VecGetSize(designDomain,&numDesign); CHKERRQ(ierr);
  MMA mma(numDesign,numConstraints,designDomain);
  PetscScalar tol = 1e-4;
  PetscScalar ch = 1e-5;
  PetscInt maxIter = 1000;
  PetscInt curIter = 0;
  PetscScalar movlim = 0.2;
  char output[50];
  Vec conSens[2];
  Vec conDesignSens[2];
  adjSolver.getSensitivities(&sens);
  presAdjSolver.getSensitivities(&presSens);
  conSens[0] = presSens;
  conSens[1] = vcSens;
  conDesignSens[0] = designPresSens;
  conDesignSens[1] = designVcSens;
  PetscScalar conValue[2];
  PetscPrintf(PETSC_COMM_WORLD,"-----------------------\n");

  // Optimization loop
  while (ch > tol && curIter < maxIter){

    ierr = map.mapDesignToFull(designDomain,fullDomain); CHKERRQ(ierr);
    ierr = filter.filterDesign(fullDomain,physicalFullDomain); CHKERRQ(ierr);
    ierr = solver.setFieldFromDesignVec(physicalFullDomain); CHKERRQ(ierr);

    if (curIter % 50 == 0){
      sprintf(output,"Design_iter_%d.vts",curIter);
      ierr = solver.outputFields(outputFolder,output); CHKERRQ(ierr);
    }

    // Compute objective and sensitivities
    ierr = opt(solver,adjSolver,presAdjSolver,objective,stateConstraint,&objValue,&presValue,
               &sens,&presSens); CHKERRQ(ierr);
    // Rescale pressure drop constraint
    presValue = presValue/(presConst*refPressure) - 1.;
    ierr = VecScale(presSens,1./(presConst*refPressure)); CHKERRQ(ierr);

    // Volume constraint
    ierr = map.mapFullToDesign(physicalFullDomain,physicalDesignDomain); CHKERRQ(ierr);
    vc.computeFunctional(physicalDesignDomain,&volumeValue);
    vc.computeSensitivities(vcSens);

    // Filter sensitivities
    ierr = filter.filterSensitivities(sens,conSens,numConstraints); CHKERRQ(ierr);
    // Map to design domain
    ierr = map.mapFullToDesign(sens,designSens); CHKERRQ(ierr);
    ierr = map.mapFullToDesign(presSens,designPresSens); CHKERRQ(ierr);
    ierr = map.mapFullToDesign(vcSens,designVcSens); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"Iter: %d, objective: %f, pressure constraint: %f, volume constraint: %f\n",
                curIter,objValue,presValue,volumeValue);
    conValue[0] = presValue;
    conValue[1] = volumeValue;
    mma.Update(designDomain,designSens,conValue,conDesignSens,movlim);
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

  // Compute final objective/constraint and output macros
  ierr = map.mapFullToDesign(physicalFullDomain,physicalDesignDomain); CHKERRQ(ierr);
  vc.computeFunctional(physicalDesignDomain,&volumeValue);
  ierr = objComp.setSteadyInit(opt.lastSteady); CHKERRQ(ierr);
  objComp(solver,objective,stateConstraint,&objValue,&presValue,true);
  presValue = presValue/(presConst*refPressure) - 1.;
  PetscPrintf(PETSC_COMM_WORLD,"Final design: objective: %f, pressure constraint: %f, volume constraint: %f\n",objValue,presValue,volumeValue);

  ierr = VecDestroy(&fullDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&physicalFullDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&vcSens); CHKERRQ(ierr);
  ierr = VecDestroy(&designDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&designSens); CHKERRQ(ierr);
  ierr = VecDestroy(&designVcSens); CHKERRQ(ierr);
  ierr = VecDestroy(&physicalDesignDomain); CHKERRQ(ierr);
  ierr = VecDestroy(&designPresSens); CHKERRQ(ierr);

  // ComputeObjective obj(25000,solver.distributionsLocal);
  // PetscScalar objValue;
  // obj(solver,objective,&objValue);
  // PetscPrintf(PETSC_COMM_WORLD,"Objective: %f\n",objValue);
  // obj(solver,objective,&objValue);
  // PetscPrintf(PETSC_COMM_WORLD,"Objective second: %f\n",objValue);

  // solver.outputFields(outputFolder,"Init.vts");
  // unsteadySolver.outputFields(outputFolder,"unsteadyInit.vts");

  // /* Main loop */
  // double t1,t2;
  // PetscInt step = 0;
  // // PetscInt numSteps = 5000;
  // Vec prev,diffVec;
  // PetscScalar procNorm;
  // PetscScalar norm = 1.;
  // char steadyOut[50];
  // ierr = VecDuplicate(solver.distributionsLocal,&prev); CHKERRQ(ierr);
  // ierr = VecDuplicate(solver.distributionsLocal,&diffVec); CHKERRQ(ierr);
  // // solver.outputMacros();
  // PetscScalar tol = 1e-6;

  // t1 = MPI_Wtime();

  // // Steady state
  // while (true){
  //   ierr = solver.streamBySwapping(); CHKERRQ(ierr);
  //   ++step;
  //   if (step % 100 == 0){
  //     ierr = VecWAXPY(diffVec,-1.,solver.distributionsLocal,prev); CHKERRQ(ierr);
  //     ierr = VecNorm(diffVec,NORM_INFINITY,&procNorm); CHKERRQ(ierr);
  //     MPI_Allreduce(&procNorm,&norm,1,MPIU_SCALAR,MPI_MAX,PETSC_COMM_WORLD);
  //     PetscPrintf(PETSC_COMM_WORLD,"Norm difference at step %d: %f\n",step,norm);
  //     // sprintf(steadyOut,"Steady_data_step_%d.vts",step);
  //     // ierr = solver.outputMacros(outputFolder,steadyOut); CHKERRQ(ierr);
  //   }
  //   if (step % 100 == 99){
  //     ierr = VecCopy(solver.distributionsLocal,prev); CHKERRQ(ierr);
  //   }
  //   ierr = solver.collideAndSwap(); CHKERRQ(ierr);

  //   if (norm < tol){
  //     break;
  //   }
  // }
  // solver.outputMacros(outputFolder,"steadyState.vts");

  // t2 = MPI_Wtime();

  // PetscPrintf(PETSC_COMM_WORLD,"Time: %f\n",t2-t1);
  // ierr = VecCopy(solver.macrosLocal,unsteadySolver.macrosLocal); CHKERRQ(ierr);

  // // Set final boundary for unsteady part
  // ierr = unsteadySolver.setDistributions(solver.distributionsLocal); CHKERRQ(ierr);
  // ierr = unsteadySolver.setCurrentTimestep(0); CHKERRQ(ierr);

  // PetscInt unsteadyStep = 0;
  // char unsteadyOut[50];
  // PetscScalar domainTemp;
  // char logname[] = "log";
  // char logpath[200] = "";
  // PetscStrcat(logpath,outputFolder);
  // PetscStrcat(logpath,logname);
  // FILE* logfile;
  // PetscFOpen(PETSC_COMM_WORLD,logpath,"w",&logfile);
  // PetscFPrintf(PETSC_COMM_WORLD,logfile,"%d %f\n",0,domainTemp);
  // PetscScalar totalTemp = 0.;
  // PetscInt numUnSteps = 25000;

  // while (unsteadyStep < numUnSteps){
  //   ierr = unsteadySolver.streamBySwapping(); CHKERRQ(ierr);
  //   totalTemp += domainTemp;
  //   ++unsteadyStep;
  //   if (unsteadyStep % 200 == 0){
  //     PetscFPrintf(PETSC_COMM_WORLD,logfile,"%d %f\n",unsteadyStep,domainTemp);
  //     // sprintf(unsteadyOut,"Unsteady_step_%d.vts",unsteadyStep);
  //     // ierr = unsteadySolver.outputMacros(outputFolder,unsteadyOut); CHKERRQ(ierr);
  //   }
  //   ierr = unsteadySolver.collideAndSwap(); CHKERRQ(ierr);
  //   unsteadySolver.computeMacroFunctional(objective,1.,&domainTemp);
  // }
  // PetscFClose(PETSC_COMM_WORLD,logfile);
  // PetscPrintf(PETSC_COMM_WORLD,"Time averaged objective: %f\n",totalTemp/numUnSteps);

}
