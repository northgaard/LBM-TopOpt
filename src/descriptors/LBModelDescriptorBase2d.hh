#ifndef LBMODELDESCRIPTORBASE2D
#define LBMODELDESCRIPTORBASE2D

#include "petsc.h"
#include "LBModelDescriptor.hh"
#include "core/geometry2d.hh"

class LBModelDescriptorBase2d : public LBModelDescriptor {

public:

  void setBoundingBoxes(const Box2d& lb, const Box2d& lbg)
  {
    localBoundingBox = lb;
    localBoundingBoxGhosted = lbg;
    hasInitializedBoxes = PETSC_TRUE;
  }

protected:

  LBModelDescriptorBase2d() : localBoundingBox(0,0,0,0),
                              localBoundingBoxGhosted(0,0,0,0), hasInitializedBoxes(PETSC_FALSE)
  {}

  Box2d localBoundingBox;
  Box2d localBoundingBoxGhosted;
  PetscBool hasInitializedBoxes;
  
};

#endif
