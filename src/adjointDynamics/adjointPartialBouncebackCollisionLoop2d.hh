#ifndef PBBADJOINTCOLLISIONLOOP2D
#define PBBADJOINTCOLLISIONLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "core/codiheader.hh"
#include "dynamics/partialBouncebackCollision.hh"
#include "adjointDynamics/reversePartialBounceback.hh"

template <class Lattice,
	  template <class> class InterpolationFunction,
	  template <class,class> class CollisionOperator,
	  template <class,class> class ReverseOperator = CollisionOperator>
class AdjointPartialBouncebackLoop2d {

public:

  AdjointPartialBouncebackLoop2d(CollisionOperator<Lattice,PetscScalar> _op,
				 ReverseOperator<Lattice,PetscScalar> _rop,
				 InterpolationFunction<PetscScalar> _inter,
				 Box2d _box)
    : boundingBox(_box), colOp(_op), reverseCol(_rop), interp(_inter) {}

  void operator()(PetscScalar*** adj, PetscScalar*** fdist, PetscScalar** obst) const
  {

    PetscScalar mac[Lattice::numMacros];
    PetscScalar macRev[Lattice::numMacros];
    PetscScalar fdistOut[Lattice::numDOF];
    PetscScalar fdistInRev[Lattice::numDOF];
    PetscScalar obstRev;
    
    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
	ReversePartialBounceback<Lattice,InterpolationFunction,
				 CollisionOperator,ReverseOperator>
	  (colOp,reverseCol,interp,fdist[jj][ii],
	   fdistOut,mac,obst[jj][ii],fdistInRev,
	   adj[jj][ii],macRev,obstRev);

	for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
	  adj[jj][ii][dd] = fdistInRev[dd];
	}
	
      }
    }
  }

private:

  Box2d boundingBox;
  CollisionOperator<Lattice,PetscScalar> colOp;
  ReverseOperator<Lattice,PetscScalar> reverseCol;
  InterpolationFunction<PetscScalar> interp;

};

/* Specialization for CoDi pack */

template <class Lattice,
	  template <class> class InterpolationFunction,
	  template <class,class> class CollisionOperator>
class AdjointPartialBouncebackLoop2d<Lattice,InterpolationFunction,
					      CollisionOperator,
					      CollisionOperator>
{

  using TapeType = ReversePetscScalar::TapeType;

public:

  AdjointPartialBouncebackLoop2d(CollisionOperator<Lattice,PetscScalar> _op,
				 InterpolationFunction<PetscScalar> _interp,
				 Box2d _box)
    : boundingBox(_box)
  {
    colOp = CollisionOperator<Lattice,PetscScalar>::
      template clone<ReversePetscScalar>(_op);
    interp = InterpolationFunction<PetscScalar>::
      template clone<ReversePetscScalar>(_interp);
  }

  void operator()(PetscScalar*** adj, PetscScalar*** fdist,
		  PetscScalar** obst, PetscScalar** os) const
  {

    ReversePetscScalar rfdistIn[Lattice::numDOF];
    ReversePetscScalar rfdistOut[Lattice::numDOF];
    ReversePetscScalar rmac[Lattice::numMacros];
    ReversePetscScalar robst;
    TapeType& tape = ReversePetscScalar::getGlobalTape();

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){

	reverseKernel(adj[jj][ii],fdist[jj][ii],obst[jj][ii],
		      rfdistIn,rfdistOut,rmac,robst,tape);

	/* Add sensitivity contribution */
	os[jj][ii] += robst.getGradient();

	/* Copy back to PETSc array */
	for (size_t id = 0; id < Lattice::numDOF; ++id){
	  adj[jj][ii][id] = rfdistIn[id].getGradient();
	}

	tape.reset();
      }
    }
  }

  void sensitivityLoop(PetscScalar*** adj, PetscScalar*** fdist,
		       PetscScalar** obst, PetscScalar** os)
  {

    ReversePetscScalar rmac[Lattice::numMacros];
    ReversePetscScalar rfdistIn[Lattice::numDOF];
    ReversePetscScalar rfdistOut[Lattice::numDOF];
    ReversePetscScalar robst;
    TapeType& tape = ReversePetscScalar::getGlobalTape();

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){

	reverseKernel(adj[jj][ii],fdist[jj][ii],obst[jj][ii],
		      rfdistIn,rfdistOut,rmac,robst,tape);

	/* Add sensitivity contribution */
	os[jj][ii] += robst.getGradient();

	tape.reset();
      }
    }
  }

private:

  Box2d boundingBox;
  CollisionOperator<Lattice,ReversePetscScalar> colOp;
  InterpolationFunction<ReversePetscScalar> interp;

  inline void reverseKernel(PetscScalar* adj, PetscScalar* fdist, PetscScalar& obst,
			    ReversePetscScalar* rfdistIn, ReversePetscScalar* rfdistOut,
			    ReversePetscScalar* rmac, ReversePetscScalar& robst,
			    TapeType& tape) const
  {

    tape.setActive();

    /* Initialize reverse values */
    for (size_t id = 0; id < Lattice::numDOF; ++id){
      rfdistIn[id] = fdist[id];
      tape.registerInput(rfdistIn[id]);
    }
    robst = obst;
    tape.registerInput(robst);

    PartialBounceback<Lattice,InterpolationFunction,
		      CollisionOperator,ReversePetscScalar>
      (colOp,interp,rfdistIn,rfdistOut,rmac,robst);

    /* Set output values */
    for (size_t id = 0; id < Lattice::numDOF; ++id){
      tape.registerOutput(rfdistOut[id]);
      rfdistOut[id].setGradient(adj[id]);
    }
	
    tape.setPassive();
    tape.evaluate();

  }

};

#endif
