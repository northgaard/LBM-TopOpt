#ifndef ADJOINTISOTHERMALLBLOOPS2D
#define ADJOINTISOTHERMALLBLOOPS2D

#include "petsc.h"
#include "core/codiheader.hh"
#include "core/geometry2d.hh"
#include "adjointDynamics/adjointLBLoop.hh"

/*
  Adjoint collision
*/

/* CoDiPack */
template <class CollisionOperator, class CodiCollisionOperator>
class CodiAdjointObstIsoCollisionLoop2d : public AdjointLBLoop {

  using Real = typename CollisionOperator::RealType;
  using Reverse = typename CodiCollisionOperator::RealType;
  using Lattice = typename CollisionOperator::LatticeType;
  using Tape = typename Reverse::TapeType;
public:
  CodiAdjointObstIsoCollisionLoop2d(Box2d _box, CodiCollisionOperator _op)
    : adjColOp(_op), boundingBox(_box){}
  void execute(void* adj_v, void* fdist_v, void* obst_v,
               void* sens_v) override
  {
    auto adj = static_cast<Real***>(adj_v);
    auto fdist =  static_cast<Real***>(fdist_v);
    auto obst = static_cast<Real***>(obst_v);
    auto sens = static_cast<Real**>(sens_v);
    Tape& tape = Reverse::getGlobalTape();

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        tape.setActive();
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          rfdistIn[dd] = fdist[jj][ii][dd];
          tape.registerInput(rfdistIn[dd]);
        }
        for (size_t dd = 0; dd < CollisionOperator::numAdditionalFields; ++dd){
          robst[dd] = obst[jj][ii][dd];
          tape.registerInput(robst[dd]);
        }

        adjColOp(rfdistIn,rfdistOut,rmac,robst);

        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          tape.registerOutput(rfdistOut[dd]);
          rfdistOut[dd].setGradient(adj[jj][ii][dd]);
        }
        tape.setPassive();
        tape.evaluate();

        sens[jj][ii] += robst[0].getGradient();

        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          adj[jj][ii][dd] = rfdistIn[dd].getGradient();
        }
        tape.reset();
      }
    }
  }
private:
  Reverse rfdistIn[Lattice::numDOF];
  Reverse rfdistOut[Lattice::numDOF];
  Reverse rmac[3];
  Reverse robst[CollisionOperator::numAdditionalFields];
  CodiCollisionOperator adjColOp;
  Box2d boundingBox;
};

/* Source transformation */
template <class CollisionOperator, class ReverseCollisionOperator>
class SourceAdjointObstIsoCollisionLoop2d : public AdjointLBLoop {

  using Real = typename CollisionOperator::RealType;
  using Lattice = typename CollisionOperator::LatticeType;
public:
  SourceAdjointObstIsoCollisionLoop2d(Box2d _box, ReverseCollisionOperator _op)
    : adjColOp(_op), boundingBox(_box){}
  void execute(void* adj_v, void* fdist_v, void* obst_v, void* sens_v) override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    auto sens = static_cast<PetscScalar**>(sens_v);
    Real fdistOut[Lattice::numDOF];
    Real rfdistIn[Lattice::numDOF];
    Real mac[3];
    Real rmac[3];
    Real robst[CollisionOperator::numAdditionalFields];

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        // Zero reverse variables
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          rfdistIn[dd] = 0.;
        }
        for (size_t dd = 0; dd < 3; ++dd){
          rmac[dd] = 0.;
        }
        for (size_t dd = 0; dd < CollisionOperator::numAdditionalFields; ++dd){
          robst[dd] = 0.;
        }
        // Reverse collision
        adjColOp(fdist[jj][ii],fdistOut,mac,obst[jj][ii],rfdistIn,
                 adj[jj][ii],rmac,robst);
        // Post processing
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          adj[jj][ii][dd] = rfdistIn[dd];
        }
        sens[jj][ii] += robst[0];
      }
    }
  }
private:
  ReverseCollisionOperator adjColOp;
  Box2d boundingBox;
};

/*
  Adjoint stream
*/

template <class CollisionOperator>
class AdjointIsoStreamingLoop2d : public AdjointLBLoop {

  using Lattice = typename CollisionOperator::LatticeType;
  using Real = typename CollisionOperator::RealType;
public:
  AdjointIsoStreamingLoop2d(Box2d _box, const Box2d& _boxGlobal)
    : boundingBox(_box)
  {
    xmin = _boxGlobal.xRange.getBeginId();
    xmax = _boxGlobal.xRange.getEndId();
    ymin = _boxGlobal.yRange.getBeginId();
    ymax = _boxGlobal.yRange.getEndId();
  }
  void execute(void* adj_v, void* adjPrev_v,void*,void*) override
  {
    auto adj = static_cast<Real***>(adj_v);
    auto adjPrev = static_cast<const Real***>(adjPrev_v);
    PetscInt dd,nextx,nexty;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        for (dd = 1; dd < Lattice::numDOF; ++dd){
          nextx = ii + Lattice::ex[dd];
          nexty = jj + Lattice::ey[dd];

          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              adj[jj][ii][dd] = adjPrev[nexty][nextx][dd];
            } else
            {
              // Default is mid-way bounceback
              adj[jj][ii][dd] = adjPrev[jj][ii][Lattice::opp[dd]];
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
