#ifndef LBSOLVERBASE
#define LBSOLVERBASE

#include "petsc.h"

class LBSolverBase {

public:

  LBSolverBase() : distributionsGlobal(nullptr), distributionsLocal(nullptr),
                   macrosGlobal(nullptr), initMacrosGlobal(nullptr), currentTimestep(0) {}
  virtual ~LBSolverBase(){}

  virtual PetscErrorCode streamAndCollide() = 0;
  virtual PetscErrorCode collideAndStream() = 0;
  virtual PetscErrorCode collide() = 0;
  virtual PetscErrorCode stream() = 0;
  virtual PetscErrorCode initializeDistributions() = 0;
  virtual PetscErrorCode outputDistributions() = 0;
  virtual PetscErrorCode outputMacros() = 0;

  PetscErrorCode setDistributions(Vec dist)
  {

    // PetscErrorCode ierr;
    // ierr = VecCopy(dist,distributionsGlobal); CHKERRQ(ierr);

    return 0;
    
  }

  void setCurrentTimestep(PetscInt _ct)
  {
    currentTimestep = _ct;
  }

  Vec distributionsGlobal;
  Vec distributionsLocal;
  Vec macrosGlobal;
  Vec initMacrosGlobal;

protected:

  PetscInt currentTimestep;

};

#endif
