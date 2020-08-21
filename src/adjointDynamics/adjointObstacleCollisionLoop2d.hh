#ifndef ADJOINTOBSTLOOP2D
#define ADJOINTOBSTLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "core/memoryHandling.hh"
#include "core/codiheader.hh"
#include "latticeBoltzmann/obstacleLBSolverBase2d.hh"

template <class CodiObstacleCollision>
class CodiAdjointObstacleCollisionLoop2d {

  template <class T>
  friend class ObstacleCollisionLoop2d;
  
  using TapeType = ReversePetscScalar::TapeType;

public:

  void operator()(PetscScalar*** adj, PetscScalar*** fdist,
		  PetscScalar** obst, PetscScalar** os) const
  {

    TapeType& tape = ReversePetscScalar::getGlobalTape();

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){

	tape.setActive();

	for (size_t id = 0; id < latticeDOF; ++id){
	  rfdistIn[id] = fdist[jj][ii][id];
	  tape.registerInput(rfdistIn[id]);
	}
	robst = obst[jj][ii];
	tape.registerInput(robst);

	codiOp(rfdistIn,rfdistOut,rmac,robst);

	for (size_t id = 0; id < latticeDOF; ++id){
	  tape.registerOutput(rfdistOut[id]);
	  rfdistOut[id].setGradient(adj[jj][ii][id]);
	}

	tape.setPassive();
	tape.evaluate();

	os[jj][ii] += robst.getGradient();

	for (size_t id = 0; id < latticeDOF; ++id){
	  adj[jj][ii][id] = rfdistIn[id].getGradient();
	}

	tape.reset();
	
      }
    }
  }

private:

  CodiAdjointObstacleCollisionLoop2d(const ObstacleLBSolverBase2d& solver,
				     CodiObstacleCollision _op)
    : boundingBox(solver.getLocalBoundingBox()), codiOp(_op)
  {

    PetscErrorCode ierr;
    ierr = DMDAGetInfo(solver.latticeGrid,0,0,0,0,0,0,0,&latticeDOF,0,0,0,0,0);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = DMDAGetInfo(solver.macroGrid,0,0,0,0,0,0,0,&macroDOF,0,0,0,0,0);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);

    ReversePetscScalar* temp;
    ierr = PetscMalloc1(2*latticeDOF + macroDOF,&temp);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = registerPointerForDeallocation(temp);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);

    // Align pointers
    rfdistIn = temp;
    rfdistOut = temp + latticeDOF;
    rmac = temp + 2*latticeDOF;
    
  }

  Box2d boundingBox;
  CodiObstacleCollision codiOp;
  mutable ReversePetscScalar* rfdistIn;
  mutable ReversePetscScalar* rfdistOut;
  mutable ReversePetscScalar* rmac;
  mutable ReversePetscScalar robst;
  PetscInt latticeDOF, macroDOF;

};

#endif
