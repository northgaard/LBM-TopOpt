#ifndef NEWLBSOLVER2D
#define NEWLBSOLVER2D

#include "petsc.h"
#include "LBSolver.hh"
#include "core/geometry2d.hh"
#include "core/definitions2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"

class NewLBSolverInfo2d {

  friend class NewLBSolver2d;

public:

  NewLBSolverInfo2d()
  {
    boundaryX = DM_BOUNDARY_NONE;
    boundaryY = DM_BOUNDARY_NONE;
  }
  PetscInt nx,ny;

private:

  DMBoundaryType boundaryX, boundaryY;

};

class NewLBSolver2d : public LBSolver {

public:

  Box2d getBoundingBox() const { return globalBoundingBox; }
  Box2d getLocalBoundingBox() const { return localBoundingBox; }
  Box2d getLocalBoundingBoxGhosted() const { return localBoundingBoxGhosted; }
  PetscErrorCode addBoundaryCondition(const Box2d,const NewBoundaryDescriptor2d&);
  PetscErrorCode addBoundaryConditionRaw(const Box2d,const NewBoundaryDescriptor2d&,
                                         BoundaryOrientation2d);
  PetscErrorCode setUniformFieldValues(Box2d,PetscScalar*);

protected:

  NewLBSolver2d(){}
  PetscErrorCode createPetscObjects(const NewLBSolverInfo2d&,const PetscInt,
                                    const PetscInt,const PetscInt);
  PetscErrorCode createBoundingBoxes(const NewLBSolverInfo2d&);

  Box2d globalBoundingBox;
  Box2d localBoundingBox;
  Box2d localBoundingBoxGhosted;

};

#endif
