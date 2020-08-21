#ifndef LBSOLVER
#define LBSOLVER

#include "petsc.h"
#include "dynamics/LBLoop.hh"
#include "boundaryConditions/LBBoundaryLoop.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"
#include "adjointLatticeBoltzmann/adjointLBSolver.hh"
#include "functionals/LBMacroFunctional.hh"
#include <functional>
#include <vector>
#include <string>

class LBSolver {

public:
  PetscErrorCode streamAndCollide();
  PetscErrorCode collideAndStream();
  PetscErrorCode collide();
  // PetscErrorCode collideAndSwap();
  PetscErrorCode stream();
  // PetscErrorCode streamBySwapping();
  PetscErrorCode computeMacros();
  PetscErrorCode initializeAtEquilibrium();
  PetscErrorCode outputDistributions(const std::string&, const std::string&);
  PetscErrorCode outputMacros(const std::string&, const std::string&);
  PetscErrorCode outputFields(const std::string&, const std::string&);
  PetscErrorCode computeMacroFunctional(const LBMacroFunctional&,
                                        const PetscScalar,PetscScalar*);
  PetscErrorCode setDistributions(Vec);
  PetscErrorCode setCurrentTimestep(PetscInt);
  PetscErrorCode getDesignGrid(DM*);
  PetscErrorCode createDesignVec(Vec*);
  PetscErrorCode setFieldFromDesignVec(Vec);

  /* Communicator */
  MPI_Comm communicator;

  /* Petsc objects */
  Vec distributionsGlobal;
  Vec macrosGlobal;
  Vec initMacrosGlobal;
  Vec materialGlobal;
  DM latticeGrid;
  DM macroGrid;
  DM materialGrid;
protected:
  LBSolver();
  ~LBSolver();
  void boundaryComputations(void*);
  PetscErrorCode initializeDesignGrid();

  LBLoop* collideAndSwapLoop;
  LBLoop* collideAndStreamLoop;
  LBLoop* collideLoop;
  LBLoop* streamBySwappingLoop;
  LBLoop* streamLoop;
  LBLoop* computeMacrosLoop;
  LBLoop* initializeAtEquilibriumLoop;
  std::vector<LBBoundaryLoop*> boundaryContainer;

  enum SolverState {STREAMED,COLLIDEDANDSWAPPED,COLLIDED,INITIALIZED,INVALID};
  SolverState currentState;
  PetscInt currentTimestep;

  /* Adjoint solver methods */
  PetscErrorCode setUpAdjointSolver(AdjointLBSolver*);
  PetscErrorCode setAdjointCollisionLoop(AdjointLBLoop*,AdjointLBSolver*);
  PetscErrorCode setAdjointStreamingLoop(AdjointLBLoop*,AdjointLBSolver*);
  PetscErrorCode addAdjointBoundaryCondition(AdjointLBBoundaryLoop*,AdjointLBSolver*);
private:
  DM designGrid;
};

#endif
