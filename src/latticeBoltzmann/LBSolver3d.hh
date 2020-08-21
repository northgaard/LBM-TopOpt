#ifndef LBSOLVER3D
#define LBSOLVER3D

#include "petsc.h"
#include "LBSolver.hh"
#include "core/geometry3d.hh"
#include "boundaryConditions/boundaryDescriptor3d.hh"

class LBSolverInfo3d {

  friend class LBSolver3d;
public:
  LBSolverInfo3d()
  {
    boundaryX = DM_BOUNDARY_NONE;
    boundaryY = DM_BOUNDARY_NONE;
    boundaryZ = DM_BOUNDARY_NONE;
  }
  PetscInt nx,ny,nz;
private:
  DMBoundaryType boundaryX, boundaryY, boundaryZ;
};

class LBSolver3d : public LBSolver {

public:
  Box3d getBoundingBox() const { return globalBoundingBox; }
  Box3d getLocalBoundingBox() const { return localBoundingBox; }
  Box3d getLocalBoundingBoxGhosted() const { return localBoundingBoxGhosted; }
  PetscErrorCode addBoundaryCondition(const Box3d, const BoundaryDescriptor3d&);
  PetscErrorCode setUniformFieldValues(Box3d,PetscScalar*);
protected:
  LBSolver3d(){}
  PetscErrorCode createPetscObjects(const LBSolverInfo3d&,const PetscInt,
                                    const PetscInt,const PetscInt);
  PetscErrorCode createBoundingBoxes(const LBSolverInfo3d&);

  Box3d globalBoundingBox;
  Box3d localBoundingBox;
  Box3d localBoundingBoxGhosted;
};

#endif
