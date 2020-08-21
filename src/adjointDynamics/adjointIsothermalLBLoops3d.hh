#ifndef ADJOINTISOTHERMALLBLOOPS3D
#define ADJOINTISOTHERMALLBLOOPS3D

#include "petsc.h"
#include "core/geometry3d.hh"
#include "adjointDynamics/adjointLBLoop.hh"
#include "traits/collisionTraits3d.hh"

/*
  Adjoint collision
*/

/* Source transformation */
template <class CollisionOperator, class ReverseCollisionOperator>
class SourceAdjointObstIsoCollisionLoop3d : public AdjointLBLoop {

  using Real = CollisionOperator::RealType;
  using Lattice = typename get_lattice<CollisionOperator>::type;
  static constexpr numAdditionalFields =
    get_num_additional_fields<CollisionOperator>::value;
public:
  SourceAdjointObstIsoCollisionLoop3d(Box3d _box, ReverseCollisionOperator _op)
    : adjColOp(_op), boundingBox(_box){}
  void execute(void* adj_v, void* fdist_v, void* obst_v, void* sens_v) override
  {
    auto adj = static_cast<PetscScalar****>(adj_v);
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto obst = static_cast<PetscScalar****>(obst_v);
    auto sens = static_cast<PetscScalar***>(sens_v);
    Real fdistOut[Lattice::numDOF];
    Real rfdistIn[Lattice::numDOF];
    Real mac[4];
    Real rmac[4];
    Real robst[numAdditionalFields];

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          // Zero reverse variables
          for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
            rfdistIn[dd] = 0.;
          }
          for (size_t dd = 0; dd < 4; ++dd){
            rmac[dd] = 0.;
          }
          for (size_t dd = 0; dd < numAdditionalFields; ++dd){
            robst[dd] = 0.;
          }
          // Reverse collision
          adjColOp(fdist[kk][jj][ii],fdistOut,mac,obst[kk][jj][ii],rfdistIn,
                   adj[kk][jj][ii],rmac,robst);
          // Post processing
          for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
            adj[kk][jj][ii][dd] = rfdistIn[dd];
          }
          sens[kk][jj][ii] += robst[0];
        }
      }
    }
  }
private:
  ReverseCollisionOperator adjColOp;
  Box3d boundingBox;
};

/*
  Adjoint stream
*/

template <class CollisionOperator>
class AdjointIsoStreamingLoop3d : public AdjointLBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
  using Real = CollisionOperator::RealType;
public:
  AdjointIsoStreamingLoop3d(Box3d _box, const Box3d& _boxGhosted)
    : boundingBox(_box)
  {
    xmin = _boxGhosted.xRange.getBeginId();
    xmax = _boxGhosted.xRange.getEndId();
    ymin = _boxGhosted.yRange.getBeginId();
    ymax = _boxGhosted.yRange.getEndId();
    zmin = _boxGhosted.zRange.getBeginId();
    zmax = _boxGhosted.zRange.getEndId();
  }
  void execute(void* adj_v, void* adjPrev_v, void*,void*) override
  {
    auto adj = static_cast<Real****>(adj_v);
    auto adjPrev = static_cast<Real****>(adjPrev_v);
    PetscInt dd,nextx,nexty,nextz;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.zRange){
          for (dd = 0; dd < Lattice::numDOF; ++dd){
            nextx = ii + Lattice::ex[dd];
            nexty = jj + Lattice::ey[dd];
            nextz = kk + Lattice::ez[dd];

            if (PetscLikely(nextx >= xmin && nexty >= ymin && nextz >= zmin &&
                            nextx <= xmax && nexty <= ymax && nextz <= zmax))
              {
                adj[kk][jj][ii][dd] = adjPrev[nextz][nexty][nextx][dd];
              } else {
              adj[kk][jj][ii][dd] = 0.;
            }
          }
        }
      }
    }
  }
private:
  Box3d boundingBox;
  PetscInt xmin,xmax;
  PetscInt ymin,ymax;
  PetscInt zmin,zmax;
};

#endif
