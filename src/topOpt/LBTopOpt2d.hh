#ifndef LBTOPOPT2D
#define LBTOPOPT2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/obstacleLBSolver2d.hh"
#include "adjointLatticeBoltzmann/adjointLBSolver2d.hh"
#include "topOpt/forwardSolutionStorage.hh"
#include "topOpt/filter.hh"
#include "topOpt/MMA.hh"
#include "constraints/volumeConstraint.hh"
#include <vector>
#include <memory>

class LBTopOpt2d {

public:
  
  LBTopOpt2d(ObstacleLBSolver2d&, AdjointLBSolver2d&,
	     std::unique_ptr<Filter> = nullptr);
  virtual ~LBTopOpt2d()
  {
    constraintValues.clear();
    constraintSensitivities.clear();
  }

  PetscErrorCode getDesignVector(Vec&);
  PetscErrorCode restoreDesignVector(const Vec&);
  PetscErrorCode getObjectiveAndConstraints(PetscScalar&,PetscScalar&);
  PetscErrorCode
  getObjectiveConstraintsAndSensitivities(PetscScalar&,PetscScalar&,
					     Vec&,std::vector<Vec>&);
  PetscErrorCode setDesignDomain(const Box2d);
  PetscErrorCode addVolumeConstraint(PetscScalar,VolumeConstraintType);
  std::unique_ptr<MMA> getMMA();

  PetscErrorCode outputFullDomain(PetscInt);
    
  Vec physicalFullDomainVec;
  Vec fullDomainVec;
  Vec designDomainVec;
  Vec physicalDesignDomainVec;
  Vec sensitivityFull;
  Vec sensitivityDesign;
  DM designGrid;
  VecScatter fullToDesign;
  ObstacleLBSolver2d& solver;
  AdjointLBSolver2d& adjointSolver;

protected:

  PetscErrorCode createDesignObjects(const Box2d&);
  PetscErrorCode forwardAllocation(PetscInt,ForwardSolutionStorage*);
  PetscErrorCode resetSensitivities();
  Vec& getGlobalStateVector(){ return solver.distributionsGlobal; }

  PetscScalar objectiveValue;
  std::vector<PetscScalar> constraintValues;
  std::vector<Vec> constraintSensitivities;
  VolumeConstraint constraint;
  std::unique_ptr<Filter> filter;

};

#endif
