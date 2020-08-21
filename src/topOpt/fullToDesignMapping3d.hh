#ifndef FULLTODESIGNMAPPING3D
#define FULLTODESIGNMAPPING3D

#include "petsc.h"
#include "core/geometry3d.hh"
#include "latticeBoltzmann/LBSolver3d.hh"

class FullToDesignMapping3d {

public:
  FullToDesignMapping3d(Box3d,LBSolver3d&);
  ~FullToDesignMapping3d();
  PetscErrorCode mapFullToDesign(Vec,Vec);
  PetscErrorCode mapDesignToFull(Vec,Vec);
  PetscErrorCode createDesignDomainVec(Vec*);
  PetscErrorCode createWorldDesignDomainVec(Vec*);
  MPI_Comm communicator;
private:
  DM designGrid;
  DM worldDesignGrid;
  VecScatter fullToDesign;
};

#endif
