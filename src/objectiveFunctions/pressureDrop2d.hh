#ifndef PRESSUREDROP2D_O
#define PRESSUREDROP2D_O

#include "petsc.h"
#include "core/geometry2d.hh"
#include "objectiveFunctions/objectiveFunction2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"
#include <memory>

template <class Lattice>
class PressureDrop2d : public ObjectiveFunction2d {

public:

  virtual void
  evaluate(PetscScalar*** mac, PetscScalar& value) const
  {
    
    value = 0.;
    PetscScalar ps = 0.;

    for (auto jj : inflow.yRange){
      for (auto ii : inflow.xRange){
        ps += invInflowLength * mac[jj][ii][0];
      }
    }

    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        ps -= invOutflowLength * mac[jj][ii][0];
      }
    }

    MPI_Allreduce(&ps,&value,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);

  }
  
  virtual void
  adjointStateSource(PetscInt,PetscScalar*** adj,
		     PetscScalar*** fdist) const
  {

    for (auto jj : inflow.yRange){
      for (auto ii : inflow.xRange){
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          adj[jj][ii][dd] += invInflowLength*invTimesteps;
        }
      }
    }

    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          adj[jj][ii][dd] -= invOutflowLength*invTimesteps;
        }
      }
    }

  }

  static std::unique_ptr<ObjectiveFunction2d>
  make (const Box2d _ginflow, Box2d _goutflow,
	const LBSolverBase2d& solver, PetscInt numTimesteps)
  {
    return std::unique_ptr<ObjectiveFunction2d>
      (new PressureDrop2d<Lattice>(_ginflow,_goutflow,
                                   solver,numTimesteps));
  }
  

private:

  PressureDrop2d(const Box2d _ginflow, const Box2d _goutflow,
		 const LBSolverBase2d& solver, PetscInt numTimesteps)
  {

    boxIntersection(solver.getLocalBoundingBox(),_ginflow,inflow);
    boxIntersection(solver.getLocalBoundingBox(),_goutflow,outflow);

    PetscInt inflowLength = _ginflow.getNx()*_ginflow.getNy();
    PetscInt outflowLength = _goutflow.getNx()*_goutflow.getNy();

    invInflowLength = (PetscScalar) (1./inflowLength);
    invOutflowLength = (PetscScalar) (1./outflowLength);
    invTimesteps = (PetscScalar) (1./numTimesteps);

  }

  Box2d inflow;
  Box2d outflow;
  PetscScalar invInflowLength,invOutflowLength,invTimesteps;
  
};

#endif
