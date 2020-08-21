#ifndef ISOLBDESCRIPTOR2D
#define ISOLBDESCRIPTOR2D

#include "petsc.h"
#include "LBModelDescriptorBase2d.hh"
#include "dynamics/collisionLoop2d.hh"
#include "dynamics/initializationLoop2d.hh"

template <class Lattice, class CollisionOperator>
class IsothermalLBDescriptor2d : public LBModelDescriptorBase2d {

public:

  IsothermalLBDescriptor2d(CollisionOperator _op) : colOp(_op){}

  PetscInt getNumDOF() override { return Lattice::numDOF; }
  PetscInt getNumMacros() override { return Lattice::numMacros; }
  PetscInt getNumAdditionalFields() override { return 0; }

  PetscErrorCode getCollideAndSwapFunc(DynamicsFunction* func) override
  {
    CollisionAndSwapLoop2d<Lattice,CollisionOperator>
      functionObject(localBoundingBoxGhosted,colOp);
    *func = functionObject;
    return 0;
  }
  PetscErrorCode getStreamBySwappingFunc(DynamicsFunction* func) override
  {
    StreamBySwappingLoop2d<Lattice> functionObject(localBoundingBoxGhosted);
    *func = functionObject;
    return 0;
  }
  PetscErrorCode getCollideAndStreamFunc(DynamicsFunction* func) override
  {
    CollideAndStreamSingleLoop2d<Lattice,CollisionOperator>
      functionObject(localBoundingBoxGhosted,colOp);
    *func = functionObject;
    return 0;
  }
  PetscErrorCode getInitializeAtEquilibriumFunc(DynamicsFunction* func) override
  {
    EquilibriumInitializationLoop2d<Equilibrium> functionObject(localBoundingBoxGhosted); 
    *func = functionObject;
    return 0;
  }

private:

  using Equilibrium = typename CollisionOperator::Equilibrium;

  CollisionOperator colOp;

};

#endif
