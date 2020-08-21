#ifndef ADJOINTTHERMALSTREAMINGLOOP2D
#define ADJOINTTHERMALSTREAMINGLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"

template <class IsoLattice, class ThermalLattice>
class AdjointThermalStreamingLoop2d {

  template <class T, class U>
  friend class ThermalStreamingLoop2d;

public:

  void operator()(PetscScalar*** lambda, PetscScalar*** lambdaPrev) const
  {

    PetscInt dd;
    PetscInt nextx,nexty;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){

	for (dd = 0; dd < IsoLattice::numDOF; ++dd){

	  nextx = ii + IsoLattice::ex[dd];
	  nexty = jj + IsoLattice::ey[dd];

	  if (PetscLikely(nextx >= xmin && nexty >= ymin &&
			  nextx <= xmax && nexty <= ymax))
	    {
	      lambda[jj][ii][dd] = lambdaPrev[nexty][nextx][dd];
	    } else {
	    lambda[jj][ii][dd] = 0.;
	  }
	}

	for (dd = 0; dd < ThermalLattice::numDOF; ++dd){

	  nextx = ii + ThermalLattice::ex[dd];
	  nexty = jj + ThermalLattice::ey[dd];

	  if (PetscLikely(nextx >= xmin && nexty >= ymin &&
			  nextx <= xmax && nexty <= ymax))
	    {
	      lambda[jj][ii][dd + IsoLattice::numDOF] =
		lambdaPrev[nexty][nextx][dd + IsoLattice::numDOF];
	    } else {
	    lambda[jj][ii][dd + IsoLattice::numDOF] = 0.;
	  }
	}
      }
    }
    
  }

private:

  AdjointThermalStreamingLoop2d(const LBSolverBase2d& solver)
    : boundingBox(solver.getLocalBoundingBox())
  {
    Box2d _boxG = solver.getLocalBoundingBoxGhosted();
    
    xmin = _boxG.xRange.getBeginId();
    xmax = _boxG.xRange.getEndId();
    ymin = _boxG.yRange.getBeginId();
    ymax = _boxG.yRange.getEndId();
  }

  Box2d boundingBox;
  PetscInt xmin,xmax;
  PetscInt ymin,ymax;

};

#endif
