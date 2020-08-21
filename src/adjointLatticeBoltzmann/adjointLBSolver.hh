#ifndef ADJOINTLBSOLVER
#define ADJOINTLBSOLVER

#include "petsc.h"
#include "adjointDynamics/adjointLBLoop.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"
#include "functionals/LBMacroFunctional.hh"
#include <vector>
#include <string>

class AdjointLBSolver {

  friend class LBSolver;
public:
  AdjointLBSolver();
  ~AdjointLBSolver();
  PetscErrorCode adjointCollide(Vec,Vec);
  PetscErrorCode adjointStream(Vec);
  PetscErrorCode computeCollideSource(const LBMacroFunctional&,const PetscScalar,
                                      Vec,Vec);
  PetscErrorCode computeSensitivitySource(const LBMacroFunctional&,const PetscScalar,
                                          Vec,Vec);
  PetscErrorCode outputAdjoints(const std::string&, const std::string&);
  PetscErrorCode resetAdjoints();
  PetscErrorCode setCurrentTimestep(PetscInt);
  PetscErrorCode getSensitivities(Vec*);

  /* MPI */
  MPI_Comm communicator;

  /* Petsc objects */
  Vec adjointGlobal;
  Vec sensitivityGlobal;
  DM latticeGrid;
  DM materialGrid;
  DM sensitivityGrid;
private:
  PetscErrorCode createPetscObjects(DM,DM,DM);
  void adjointBoundaryComputations(void*,void*);

  AdjointLBLoop* adjointCollideLoop;
  AdjointLBLoop* adjointStreamLoop;
  std::vector<AdjointLBBoundaryLoop*> adjointBoundaryContainer;

  PetscInt currentTimestep;
};

#endif
