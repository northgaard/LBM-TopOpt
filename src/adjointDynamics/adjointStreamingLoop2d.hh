#ifndef ADJSTREAMINGLOOP2D
#define ADJSTREAMINGLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"

template <class Lattice>
class AdjointStreamingLoop2d {

public:

  AdjointStreamingLoop2d(Box2d,const Box2d&);
  void operator()(PetscScalar***,PetscScalar***) const;

private:

  Box2d boundingBox;
  PetscInt xmin,xmax;
  PetscInt ymin,ymax;

};

template <class Lattice>
AdjointStreamingLoop2d<Lattice>::AdjointStreamingLoop2d(Box2d _box,
							const Box2d& _boxG)
  : boundingBox(_box)
{
  xmin = _boxG.xRange.getBeginId();
  xmax = _boxG.xRange.getEndId();
  ymin = _boxG.yRange.getBeginId();
  ymax = _boxG.yRange.getEndId();
}

template <class Lattice>
void AdjointStreamingLoop2d<Lattice>::operator()(PetscScalar*** lambda,
						 PetscScalar*** lambdaPrev) const
{

  PetscInt dd;
  PetscInt nextx,nexty;

  for (auto jj : boundingBox.yRange){
    for (auto ii : boundingBox.xRange){
      for (dd = 0; dd < Lattice::numDOF; ++dd){

	// Streaming is opposite direction of forward solution
	nextx = ii + Lattice::ex[dd];
	nexty = jj + Lattice::ey[dd];

	if (PetscLikely(nextx >= xmin && nexty >= ymin &&
			nextx <= xmax && nexty <= ymax))
	  {
	    lambda[jj][ii][dd] = lambdaPrev[nexty][nextx][dd];
	  } else {
	  lambda[jj][ii][dd] = 0.;
	}
      }
    }
  }

}

#endif
