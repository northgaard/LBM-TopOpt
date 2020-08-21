#ifndef BASECOLLISIONLOOP2D
#define BASECOLLISIONLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"

template <class CollisionOperator>
class BaseCollisionLoop2d {

public:

  BaseCollisionLoop2d(const LBSolverBase2d& solver , CollisionOperator _op)
    : boundingBox(solver.getLocalBoundingBox()), colOp(_op) {}
  
  void operator()(PetscScalar*** fdist, PetscScalar*** mac) const
  {
    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
	colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii]); 
      }
    }
  }

private:

  Box2d boundingBox;
  CollisionOperator colOp;

};

#endif
