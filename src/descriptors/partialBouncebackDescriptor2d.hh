#ifndef PARTIALBOUNCEBACKDESCRIPTOR2D
#define PARTIALBOUNCEBACKDESCRIPTOR2D

#include "petsc.h"
#include "descriptors/collisionDescriptorBases2d.hh"
#include "dynamics/initializationLoop2d.hh"
#include "dynamics/streamingLoop2d.hh"
#include "dynamics/partialBouncebackCollisionLoop2d.hh"
#include "adjointDynamics/adjointStreamingLoop2d.hh"
#include "adjointDynamics/adjointPartialBouncebackCollisionLoop2d.hh"
#include <memory>

template <class Lattice,
	  template <class> class InterpolationFunction,
	  template <class,class> class CollisionOperator>
class PartialBouncebackDescriptor2d : public ObstacleCollisionDescriptorBase2d {

  using Operator = CollisionOperator<Lattice,PetscScalar>;
  using Equilibrium = typename Operator::Equilibrium;

public:

  virtual PetscErrorCode getInitializationLoop(Box2d localBox,
					       BaseDynamicsFunction& func) const
  {

    EquilibriumInitializationLoop2d<Equilibrium> init(localBox);
    func = [init](PetscScalar*** fdist, PetscScalar*** mac)
      { init(fdist,mac); };
    return 0;

  }

  virtual PetscErrorCode getStreamingLoop(Box2d localBox, Box2d localBoxGhosted,
					  BaseDynamicsFunction& func) const
  {

    StreamingLoop2d<Lattice> strLoop(localBox,localBoxGhosted);
    func = [strLoop](PetscScalar*** fdist, PetscScalar*** fcol)
      { strLoop(fdist,fcol); };
    return 0;

  }

  virtual PetscErrorCode getAdjointStreamingLoop(Box2d localBox,
						 Box2d localBoxGhosted,
						 BaseDynamicsFunction& func) const
  {

    AdjointStreamingLoop2d<Lattice> adjStrLoop(localBox,localBoxGhosted);
    func = [adjStrLoop](PetscScalar*** adj, PetscScalar*** adjPrev)
      { adjStrLoop(adj,adjPrev); };
    return 0;

  }

  virtual PetscErrorCode getCollisionLoop(Box2d theBox,
					  ObstacleDynamicsFunction& func) const
  {

    PartialBouncebackCollisionLoop2d<Lattice,InterpolationFunction,
				     CollisionOperator>
      colLoop(theBox,colOp,interp);
    func = [colLoop](PetscScalar*** fdist, PetscScalar*** mac, PetscScalar** obst)
      { colLoop(fdist,mac,obst); };
    return 0;
    
  }

  virtual PetscErrorCode
  getAdjointCollisionLoop(Box2d theBox,
			  AdjointObstacleDynamicsFunction& func) const
  {

    AdjointPartialBouncebackLoop2d<Lattice,InterpolationFunction,
				   CollisionOperator>
      adjColLoop(colOp,interp,theBox);
    func = [adjColLoop](PetscScalar*** adj, PetscScalar*** fdist,
			PetscScalar** obst, PetscScalar** os)
      { adjColLoop(adj,fdist,obst,os); };
    return 0;
    
  }

  virtual PetscInt getDistributionDOF() const { return Lattice::numDOF; }
  virtual PetscInt getMacroDOF() const { return Lattice::numMacros; }

  static std::unique_ptr<ObstacleCollisionDescriptorBase2d>
  make(Operator _col, InterpolationFunction<PetscScalar> _inter)
  {
    return std::unique_ptr<ObstacleCollisionDescriptorBase2d>
      (new PartialBouncebackDescriptor2d<Lattice,InterpolationFunction,
       CollisionOperator>(_col,_inter));
  }

private:

  PartialBouncebackDescriptor2d(Operator _col,
				InterpolationFunction<PetscScalar> _inter)
    : colOp(_col), interp(_inter) {}
  Operator colOp;
  InterpolationFunction<PetscScalar> interp;

};

#endif
