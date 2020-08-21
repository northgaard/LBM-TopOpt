#include "unsteadyLBTopOpt2d.hh"

PetscErrorCode UnsteadyLBTopOpt2d::computeObjective()
{

  PetscErrorCode ierr;

  /* Filter design */
  if (filter){
    ierr = filter->filterDesign(fullDomainVec,physicalFullDomainVec);
    CHKERRQ(ierr);
  } else {
    ierr = VecCopy(fullDomainVec,physicalFullDomainVec); CHKERRQ(ierr);
  }
  ierr = solver.setDomain(physicalFullDomainVec); CHKERRQ(ierr);

  ierr = solver.initializeDistributions(); CHKERRQ(ierr);
  objectiveValue = 0.;

  PetscScalar tvalue;

  for (PetscInt tt = 0; tt < numTimesteps; ++tt){
    ierr = solver.streamAndCollide(); CHKERRQ(ierr);
    ierr = solver.computeMacroObjective(obj,tvalue); CHKERRQ(ierr);
    objectiveValue += tvalue;
  }

  objectiveValue /= numTimesteps;
  
  return 0;

}

PetscErrorCode UnsteadyLBTopOpt2d::computeObjectiveWithCheckpoints()
{

  PetscErrorCode ierr;

  /* Filter design */
  if (filter){
    ierr = filter->filterDesign(fullDomainVec,physicalFullDomainVec);
    CHKERRQ(ierr);
  } else {
    ierr = VecCopy(fullDomainVec,physicalFullDomainVec); CHKERRQ(ierr);
  }
  ierr = solver.setDomain(physicalFullDomainVec); CHKERRQ(ierr);

  ierr = solver.initializeDistributions(); CHKERRQ(ierr);
  objectiveValue = 0.;
  PetscInt timestep = 0;

  // Save initial state
  ierr = store->saveForwardSolution(timestep,getGlobalStateVector()); CHKERRQ(ierr);

  PetscScalar tvalue;

  while (timestep < numTimesteps){
    ierr = solver.stream(); CHKERRQ(ierr);
    ++timestep;
    
    // Save current step
    ierr = store->saveForwardSolution(timestep,getGlobalStateVector());
    CHKERRQ(ierr);
    
    ierr = solver.collide(); CHKERRQ(ierr);
    ierr = solver.computeMacroObjective(obj,tvalue); CHKERRQ(ierr);
    objectiveValue += tvalue;

  }

  objectiveValue /= numTimesteps;
  
  return 0;
  
}

PetscErrorCode UnsteadyLBTopOpt2d::forwardAndAdjoint()
{

  PetscErrorCode ierr;

  ierr = computeObjectiveWithCheckpoints(); CHKERRQ(ierr);
  ierr = adjointSolver.resetAdjoints(); CHKERRQ(ierr);
  ierr = resetSensitivities(); CHKERRQ(ierr);

  PetscInt timestep = numTimesteps;
  Vec forwardStore;

  while (timestep > 0){

    store->getForwardSolution(timestep,forwardStore);
    ierr = adjointSolver.adjointCollide(forwardStore,solver.obstacleGlobal,
					sensitivityFull,obj.get());
    CHKERRQ(ierr);
    ierr = adjointSolver.adjointStream(); CHKERRQ(ierr);
    --timestep;

  }

  if (filter){
    ierr = filter->filterSensitivities(sensitivityFull); CHKERRQ(ierr);
  }

  return 0;

}

PetscErrorCode UnsteadyLBTopOpt2d::
allocateMemoryForForwardSolution(PetscInt numCheckPoints)
{

  PetscErrorCode ierr;

  if (numCheckPoints == -1){
    ForwardSolutionStorage* temp;
    ierr = PetscMalloc1(1,&temp); CHKERRQ(ierr);
    ierr = registerPointerForDeallocation(temp); CHKERRQ(ierr);
    ierr = forwardAllocation(numTimesteps+1,temp);
    store = temp;
  }

  return 0;

}

PetscErrorCode UnsteadyLBTopOpt2d::deallocateMemoryForForwardSolution()
{

  PetscErrorCode ierr;
  ierr = store->deallocate(); CHKERRQ(ierr);

  return 0;

}
