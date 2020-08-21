#ifndef ISOTHERMALCOLLISIONDESCRIPTOR2D
#define ISOTHERMALCOLLISIONDESCRIPTOR2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "descriptors/collisionDescriptorBases2d.hh"
#include "dynamics/initializationLoop2d.hh"
#include "dynamics/streamingLoop2d.hh"
#include "dynamics/collisionLoop2d.hh"
#include <memory>

template <class Lattice, template <class,class> class CollisionOperator>
class CollisionDescriptor2d : public CollisionDescriptorBase2d {

  using Operator = CollisionOperator<Lattice,PetscScalar>;
  using Equilibrium = typename Operator::Equilibrium;

public:

  virtual PetscErrorCode getInitializationLoop(Box2d theBox,
					       BaseDynamicsFunction& func) const
  {

    EquilibriumInitializationLoop2d<Equilibrium> init(theBox);
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

  virtual PetscErrorCode getCollisionLoop(Box2d theBox,
					  BaseDynamicsFunction& func) const
  {

    Operator colOp(par);
    CollisionLoop2d<Operator> colLoop(theBox,colOp);
    func = [colLoop](PetscScalar*** fdist, PetscScalar*** mac)
      { colLoop(fdist,mac); };
    return 0;

  }

  virtual PetscInt getDistributionDOF() const { return Lattice::numDOF; }
  virtual PetscInt getMacroDOF() const { return Lattice::numMacros; }

  static std::unique_ptr<CollisionDescriptorBase2d>
  make(IncompressibleFlowParameters par)
  {
    return std::unique_ptr<CollisionDescriptorBase2d>
      (new CollisionDescriptor2d<Lattice,CollisionOperator>(par));
  }

private:

  CollisionDescriptor2d(IncompressibleFlowParameters _p) : par(_p) {}
  IncompressibleFlowParameters par;

};

#endif
