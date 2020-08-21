#ifndef OBSTACLECOLLISIONLOOP2D
#define OBSTACLECOLLISIONLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "core/meta.hh"
#include "latticeBoltzmann/obstacleLBSolverBase2d.hh"
#include "adjointDynamics/adjointObstacleCollisionLoop2d.hh"
#include <functional>
#include <type_traits>

template <class ObstacleCollisionOperator>
class ObstacleCollisionLoop2d {

public:

  ObstacleCollisionLoop2d(const ObstacleLBSolverBase2d& solver,
			  ObstacleCollisionOperator _op)
    : boundingBox(solver.getLocalBoundingBox()), colOp(_op) {}

  void operator()(PetscScalar*** fdist, PetscScalar*** mac,
		  PetscScalar** obst) const
  {
    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
	colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
      }
    }
  }

  using AdjointDynamicsFunction = std::function<void(PetscScalar***,
						     PetscScalar***,
						     PetscScalar**,
						     PetscScalar**)>;

  template <class T = ObstacleCollisionOperator>
  AdjointDynamicsFunction getAdjoint(const ObstacleLBSolverBase2d& solver)
  {

    static_assert(has_reverse_type<ObstacleCollisionOperator>::value
		  || has_codireverse_type<ObstacleCollisionOperator>::value,
		  "The provided obstacle collision operator does not provide typedefs to any reverse operator.\n");

    if (has_reverse_type<ObstacleCollisionOperator>::value){
      return getSourceAdjoint(solver);
    } else if (has_codireverse_type<ObstacleCollisionOperator>::value){
      return getCodiAdjoint(solver);
    }

    return 0;
    
  }

private:

  template <class T = ObstacleCollisionOperator>
  AdjointDynamicsFunction getCodiAdjoint(const ObstacleLBSolverBase2d& solver,
					 typename std::enable_if<
					 has_codireverse_type<T>::value>::
					 type* = 0)
  {
    static_assert(!std::is_same<typename T::CodiReverseType,
		  void>::value,
		  "The obstacle collision operator was not initialized with the proper floating point type (PetscScalar).\n");

    using Reverse = typename T::CodiReverseType;

    Reverse reverseOp(colOp);
    CodiAdjointObstacleCollisionLoop2d<Reverse> adjointLoop(solver,reverseOp);

    return adjointLoop;
  }
  template <class T = ObstacleCollisionOperator>
  AdjointDynamicsFunction getCodiAdjoint(const ObstacleLBSolverBase2d&,
					 T* = 0, typename std::enable_if<
					 !(has_codireverse_type<T>::value)>::type* = 0)
  {
    return 0;
  }

  template <class T = ObstacleCollisionOperator>
  AdjointDynamicsFunction getSourceAdjoint(const ObstacleLBSolverBase2d& solver,
					   typename std::enable_if<
					   has_reverse_type<T>::value>::
					   type* = 0)
  {
    return 0;
  }
  template <class T = ObstacleCollisionOperator>
  AdjointDynamicsFunction getSourceAdjoint(const ObstacleLBSolverBase2d&,
					   T* = 0, typename std::enable_if<
					   !(has_reverse_type<T>::value)>::type* = 0)
  {
    return 0;
  }

  Box2d boundingBox;
  ObstacleCollisionOperator colOp;

};

#endif
