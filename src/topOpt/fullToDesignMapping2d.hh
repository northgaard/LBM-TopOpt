#ifndef FULLTODESIGNMAPPING2D
#define FULLTODESIGNMAPPING2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

class FullToDesignMapping2d {

public:
  FullToDesignMapping2d(Box2d,NewLBSolver2d&);
  ~FullToDesignMapping2d();
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
