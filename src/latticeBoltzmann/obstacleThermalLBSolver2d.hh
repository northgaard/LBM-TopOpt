#ifndef OBSTACLETHERMALLBSOLVER2D
#define OBSTACLETHERMALLBSOLVER2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "core/definitions2d.hh"
#include "dynamics/thermalStreamingLoop2d.hh"
#include "dynamics/thermalInitializationLoop2d.hh"
#include "dynamics/obstacleCollisionLoop2d.hh"
#include "latticeBoltzmann/obstacleLBSolverBase2d.hh"

class ObstacleThermalLBSolver2d : public ObstacleLBSolverBase2d {

public:

  ObstacleThermalLBSolver2d() : ObstacleLBSolverBase2d() {}
  virtual ~ObstacleThermalLBSolver2d(){}

  PetscErrorCode uniformMacroInitialization(const ThermalMacros2d&);
  PetscErrorCode getInitialMacroArray(const Box2d&,Box2d&,ThermalMacros2d***);
  PetscErrorCode restoreInitialMacroArray(const Box2d&,Box2d&,ThermalMacros2d***);

  virtual PetscErrorCode streamAndCollide();
  virtual PetscErrorCode collide();

  template <class IsoLattice, class ThermalLattice,
	    class ObstacleOperator>
  static PetscErrorCode initializeSolver(PetscInt nx, PetscInt ny,
					 ObstacleOperator op,
					 ObstacleThermalLBSolver2d& solver)
  {

    PetscErrorCode ierr;
    constexpr PetscInt numDOF = IsoLattice::numDOF + ThermalLattice::numDOF;
    constexpr PetscInt numMacros = 4;

    ierr = solver.createLatticeObjects(nx,ny,numDOF); CHKERRQ(ierr);
    ierr = solver.createMacroObjects(numMacros); CHKERRQ(ierr);
    ierr = solver.createObstacleObjects(); CHKERRQ(ierr);
    ierr = solver.createBoundingBoxes(nx,ny);

    ierr = DMDASetFieldName(solver.macroGrid,0,"rho"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(solver.macroGrid,1,"ux"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(solver.macroGrid,2,"uy"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(solver.macroGrid,3,"T"); CHKERRQ(ierr);

    // Set streaming
    ThermalStreamingLoop2d<IsoLattice,ThermalLattice>
      strLoop(solver.localBoundingBox,solver.localBoundingBoxGhosted);
    solver.streamFunc = [strLoop](PetscScalar*** fdist, PetscScalar*** fcol)
      { strLoop(fdist,fcol); };

    // Set collision
    ObstacleCollisionLoop2d<ObstacleOperator>
      opLoop(solver.localBoundingBox,op);
    solver.collisionFunc = [opLoop](PetscScalar*** fdist, PetscScalar*** mac,
				    PetscScalar** obst) 
      { opLoop(fdist,mac,obst); };

    return 0;

  }

  template <class IsoLattice,
	    class IsoEquilibrium, class ThermalEquilibrium>
  static PetscErrorCode setInitializationLoop(ObstacleThermalLBSolver2d& solver)
  {
    
    ThermalEquilibriumInitializationLoop2d
      <IsoLattice,IsoEquilibrium,ThermalEquilibrium>
      initLoop(solver.localBoundingBox);
    solver.initialization = [initLoop](PetscScalar*** fdist, PetscScalar*** mac)
      { initLoop(fdist,mac); };

    return 0;
    
  }

private:

  ObstacleDynamicsFunction collisionFunc;

};

#endif
