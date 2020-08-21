#ifndef STREAMINGLOOP2D
#define STREAMINGLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"

template <class Lattice>
class StreamingLoop2d {

public:

  StreamingLoop2d(const LBSolverBase2d& solver)
    : boundingBox(solver.getLocalBoundingBox())
  {
    Box2d _boxG = solver.getLocalBoundingBoxGhosted();

    xmin = _boxG.xRange.getBeginId();
    xmax = _boxG.xRange.getEndId();
    ymin = _boxG.yRange.getBeginId();
    ymax = _boxG.yRange.getEndId();
  }
  
  void operator()(PetscScalar*** fdist,PetscScalar*** fcol) const
  {
    PetscInt dd;
    PetscInt nextx,nexty;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        for (dd = 0; dd < Lattice::numDOF; ++dd){

          nextx = ii - Lattice::ex[dd];
          nexty = jj - Lattice::ey[dd];

          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              fdist[jj][ii][dd] = fcol[nexty][nextx][dd];
            }
        }
      }
    }
  }
  
private:

  Box2d boundingBox;
  PetscInt xmin,xmax;
  PetscInt ymin,ymax;

};

#endif
