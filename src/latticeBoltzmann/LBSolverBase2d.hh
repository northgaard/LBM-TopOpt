#ifndef LBSOLVERBASE2D
#define LBSOLVERBASE2D

#include "petsc.h"
#include "latticeBoltzmann/LBSolverBase.hh"
#include "core/geometry2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "objectiveFunctions/objectiveFunction2d.hh"
#include <functional>
#include <vector>
#include <memory>
#include <type_traits>

#ifdef PETSC35
using PetscDMBoundary = DMBoundaryType;
#else
using PetscDMBoundary = DMDABoundaryType;
#endif

class LBSolverBase2d : public LBSolverBase {

protected:

  using BaseDynamicsFunction = std::function<void(PetscScalar***,PetscScalar***)>;
  using BoundaryFunction = std::function<void(PetscInt,PetscScalar***)>;

  BaseDynamicsFunction streamFunc;
  BaseDynamicsFunction initialization;
  std::vector<BoundaryFunction> boundaryContainer;

public:

  DM latticeGrid;
  DM macroGrid;

protected:

  Box2d globalBoundingBox;
  Box2d localBoundingBox;
  Box2d localBoundingBoxGhosted;

public:

  /* Functionality */

  PetscErrorCode stream() override;
  PetscErrorCode initializeDistributions() override;
  PetscErrorCode computeMacroObjective(const std::unique_ptr<ObjectiveFunction2d>&,
				       PetscScalar&);
  
  PetscErrorCode outputMacros() override; 
  PetscErrorCode outputDistributions() override; 

  const Box2d& getBoundingBox() const { return globalBoundingBox; }
  const Box2d& getLocalBoundingBox() const { return localBoundingBox; }
  const Box2d& getLocalBoundingBoxGhosted() const
  {
    return localBoundingBoxGhosted;
  }

  /* Set up methods */
  
  inline void setStreamingFunction(BaseDynamicsFunction _fun)
  {
    streamFunc = _fun;
  }

  inline void setInitializationFunction(BaseDynamicsFunction _fun)
  {
    initialization = _fun;
  }
  
  PetscErrorCode addBoundaryCondition(const Box2d,
				      const std::unique_ptr<BoundaryDescriptor2d>&);

  template <class MacroType>
  PetscErrorCode uniformMacroInitialization(const MacroType&);
  template <class MacroType>
  PetscErrorCode getInitialMacroArray(const Box2d&,Box2d&,
				      MacroType***);
  template <class MacroType>
  PetscErrorCode restoreInitialMacroArray(const Box2d&,Box2d&,
					  MacroType***);
  
protected:

  LBSolverBase2d() : LBSolverBase(),
		     latticeGrid(0), macroGrid(0) {}
  ~LBSolverBase2d(){}

  void boundaryComputations(PetscInt,PetscScalar***);
  PetscErrorCode createLatticeObjects(PetscInt,PetscInt,PetscInt,
				      PetscDMBoundary,PetscDMBoundary);
  PetscErrorCode createMacroObjects(PetscInt);
  PetscErrorCode createBoundingBoxes(PetscInt,PetscInt);

};

/* Template implementations */

template <class MacroType>
PetscErrorCode LBSolverBase2d::
uniformMacroInitialization(const MacroType& values)
{

  static_assert(std::is_pod<MacroType>::value,"Macro type passed to uniformMacroInitialization must be POD (plain old data).\n");

  PetscErrorCode ierr;
  PetscInt dof,structSize,dd;
  const PetscScalar* valuesPointer = reinterpret_cast<const PetscScalar*>(&values);
  PetscScalar*** initMac;

  PetscFunctionBeginUser;

  structSize = sizeof(MacroType) / sizeof(PetscScalar);
  ierr = DMDAGetInfo(macroGrid,0,0,0,0,0,0,0,&dof,0,0,0,0,0); CHKERRQ(ierr);
  if (dof != structSize){
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,
	     "The size of the macro struct %d does not match the number of DOF allocated in the macro DMDA %d.\n",structSize,dof);
  }

  ierr = DMDAVecGetArrayDOF(macroGrid,initMacrosGlobal,&initMac); CHKERRQ(ierr);

  for (auto jj : localBoundingBox.yRange){
    for (auto ii : localBoundingBox.xRange){
      for (dd = 0; dd < structSize; ++dd){
        initMac[jj][ii][dd] = valuesPointer[dd];
      }
    }
  }

  ierr = DMDAVecRestoreArrayDOF(macroGrid,initMacrosGlobal,&initMac); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

template <class MacroType>
PetscErrorCode LBSolverBase2d::
getInitialMacroArray(const Box2d& inputBox,
		     Box2d& outLocalBox,
		     MacroType*** theArray)
{

  static_assert(std::is_pod<MacroType>::value,"Macro type passed to getInitialMacroArray must be POD (plain old data).\n");

  PetscErrorCode ierr;
  PetscInt dof,structSize;

  PetscFunctionBeginUser;

  structSize = sizeof(MacroType) / sizeof(PetscScalar);
  ierr = DMDAGetInfo(macroGrid,0,0,0,0,0,0,0,&dof,0,0,0,0,0); CHKERRQ(ierr);
  if (dof != structSize){
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,
	     "The size of the macro struct %d does not match the number of DOF allocated in the macro DMDA %d.\n",structSize,dof);
  }

  boxIntersection(inputBox,localBoundingBox,outLocalBox);
  ierr = DMDAVecGetArray(macroGrid,initMacrosGlobal,theArray); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

template <class MacroType>
PetscErrorCode LBSolverBase2d::
restoreInitialMacroArray(const Box2d& inputBox,
			 Box2d& outLocalBox,
			 MacroType*** theArray)
{

  static_assert(std::is_pod<MacroType>::value,"Macro type passed to restoreInitialMacroArray must be POD (plain old data).\n");

  PetscErrorCode ierr;
  PetscInt dof,structSize;

  PetscFunctionBeginUser;

  structSize = sizeof(MacroType) / sizeof(PetscScalar);
  ierr = DMDAGetInfo(macroGrid,0,0,0,0,0,0,0,&dof,0,0,0,0,0); CHKERRQ(ierr);
  if (dof != structSize){
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,
	     "The size of the macro struct %d does not match the number of DOF allocated in the macro DMDA %d.\n",structSize,dof);
  }

  ierr = DMDAVecRestoreArray(macroGrid,initMacrosGlobal,theArray); CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

#endif
