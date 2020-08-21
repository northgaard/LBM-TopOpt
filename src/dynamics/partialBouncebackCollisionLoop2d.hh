#ifndef PBBCOLLISIONLOOP2D
#define PBBCOLLISIONLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "dynamics/partialBouncebackCollision.hh"

template <class Lattice,
	  template <class> class InterpolationFunction,
	  template <class,class> class CollisionOperator>
class PartialBouncebackCollisionLoop2d {

public:

  PartialBouncebackCollisionLoop2d(Box2d _box,
				   CollisionOperator<Lattice,PetscScalar> _op,
				   InterpolationFunction<PetscScalar> _in)
    : boundingBox(_box), colOp(_op), interp(_in){}

  void operator()(PetscScalar*** fdist, PetscScalar*** mac, PetscScalar** obst) const
  {
    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
	PartialBounceback<Lattice,InterpolationFunction,
			  CollisionOperator,PetscScalar>
	  (colOp,interp,fdist[jj][ii],fdist[jj][ii],
	   mac[jj][ii],obst[jj][ii]);
      }
    }
  }

private:

  Box2d boundingBox;
  CollisionOperator<Lattice,PetscScalar> colOp;
  InterpolationFunction<PetscScalar> interp;

};

#endif
