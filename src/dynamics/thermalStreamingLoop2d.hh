#ifndef THERMALSTREAMINGLOOP2D
#define THERMALSTREAMINGLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"
#include "adjointDynamics/adjointThermalStreamingLoop2d.hh"
#include <functional>

template <class IsoLattice, class ThermalLattice>
class ThermalStreamingLoop2d {

public:

  ThermalStreamingLoop2d(const LBSolverBase2d& solver)
    : boundingBox(solver.getLocalBoundingBox())
  {

    Box2d _boxG = solver.getLocalBoundingBoxGhosted();
    
    xmin = _boxG.xRange.getBeginId();
    xmax = _boxG.xRange.getEndId();
    ymin = _boxG.yRange.getBeginId();
    ymax = _boxG.yRange.getEndId();
  }

  void operator()(PetscScalar*** fdist, PetscScalar*** fcol) const
  {

    PetscInt dd;
    PetscInt nextx,nexty;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){

	for (dd = 0; dd < IsoLattice::numDOF; ++dd){

	  nextx = ii - IsoLattice::ex[dd];
	  nexty = jj - IsoLattice::ey[dd];

	  if (PetscLikely(nextx >= xmin && nexty >= ymin &&
			  nextx <= xmax && nexty <= ymax))
	    {
	      fdist[jj][ii][dd] = fcol[nexty][nextx][dd];
	    }
	}

	for (dd = 0; dd < ThermalLattice::numDOF; ++dd){

	  nextx = ii - ThermalLattice::ex[dd];
	  nexty = jj - ThermalLattice::ey[dd];

	  if (PetscLikely(nextx >= xmin && nexty >= ymin &&
			  nextx <= xmax && nexty <= ymax))
	    {
	      fdist[jj][ii][dd + IsoLattice::numDOF] =
		fcol[nexty][nextx][dd + IsoLattice::numDOF];
	    }
	}
      }
    }

  }

  using DynamicsFunction = std::function<void(PetscScalar***,
					      PetscScalar***)>;

  // Templating so that it only compiles if method is called
  template <class T = void>
  DynamicsFunction getAdjoint(const LBSolverBase2d& solver)
  {
    AdjointThermalStreamingLoop2d<IsoLattice,ThermalLattice>
      adjLoop(solver);
    return adjLoop;
  }

private:

  Box2d boundingBox;
  PetscInt xmin,xmax;
  PetscInt ymin,ymax;

};

#endif
