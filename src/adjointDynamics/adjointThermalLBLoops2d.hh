#ifndef ADJOINTTHERMALLBLOOPS2D
#define ADJOINTTHERMALLBLOOPS2D

#include "petsc.h"
#include "core/codiheader.hh"
#include "core/geometry2d.hh"
#include "adjointDynamics/adjointLBLoop.hh"

/*
  Adjoint collision
*/

/* CoDiPack */
template <class CollisionOperator, class CodiCollisionOperator>
class CodiAdjointObstThermalCollisionLoop2d : public AdjointLBLoop {

  using Real = typename CollisionOperator::RealType;
  using Reverse = typename CodiCollisionOperator::RealType;
  using IsoLattice = typename CollisionOperator::IsoLatticeType;
  using ThermalLattice = typename CollisionOperator::ThermalLatticeType;
  using Tape = typename Reverse::TapeType;
  static constexpr PetscInt dof = IsoLattice::numDOF + ThermalLattice::numDOF;
public:
  CodiAdjointObstThermalCollisionLoop2d(Box2d _box, CodiCollisionOperator _op)
    : adjColOp(_op), boundingBox(_box){}
  void execute(void* adj_v, void* fdist_v, void* obst_v,
               void* sens_v) override
  {
    auto adj = static_cast<Real***>(adj_v);
    auto fdist = static_cast<Real***>(fdist_v);
    auto obst = static_cast<Real***>(obst_v);
    auto sens = static_cast<Real**>(sens_v);
    Tape& tape = Reverse::getGlobalTape();
    Reverse rfdistIn[dof];
    Reverse rfdistOut[dof];
    Reverse rmac[4];
    Reverse robst[CollisionOperator::numAdditionalFields];

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        tape.setActive();
        for (size_t dd = 0; dd < dof; ++dd){
          rfdistIn[dd] = fdist[jj][ii][dd];
          tape.registerInput(rfdistIn[dd]);
        }
        for (size_t dd = 0; dd < CollisionOperator::numAdditionalFields; ++dd){
          robst[dd] = obst[jj][ii][dd];
          tape.registerInput(robst[dd]);
        }

        adjColOp(rfdistIn,rfdistOut,rmac,robst);

        for (size_t dd = 0; dd < dof; ++dd){
          tape.registerOutput(rfdistOut[dd]);
          rfdistOut[dd].setGradient(adj[jj][ii][dd]);
        }
        tape.setPassive();
        tape.evaluate();

        sens[jj][ii] += robst[0].getGradient();

        for (size_t dd = 0; dd < dof; ++dd){
          adj[jj][ii][dd] = rfdistIn[dd].getGradient();
        }
        tape.reset();
      }
    }
  }
private:
  CodiCollisionOperator adjColOp;
  Box2d boundingBox;
};

template <class CollisionOperator, class ReverseCollisionOperator>
class SourceAdjointObstThermalCollisionLoop2d : public AdjointLBLoop {

  using Real = typename CollisionOperator::RealType;
  using IsoLattice = typename CollisionOperator::IsoLatticeType;
  using ThermalLattice = typename CollisionOperator::ThermalLatticeType;
  static constexpr PetscInt dof = IsoLattice::numDOF + ThermalLattice::numDOF;
public:
  SourceAdjointObstThermalCollisionLoop2d(Box2d _box, ReverseCollisionOperator _op)
    : adjColOp(_op), boundingBox(_box){}
  void execute(void* adj_v, void* fdist_v, void* obst_v, void* sens_v) override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    auto sens = static_cast<PetscScalar**>(sens_v);
    Real fdistOut[dof];
    Real rfdistIn[dof];
    Real mac[4];
    Real rmac[4];
    Real robst[CollisionOperator::numAdditionalFields];

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        // Zero reverse variables
        for (size_t dd = 0; dd < dof; ++dd){
          rfdistIn[dd] = 0.;
        }
        for (size_t dd = 0; dd < 4; ++dd){
          rmac[dd] = 0.;
        }
        for (size_t dd = 0; dd < CollisionOperator::numAdditionalFields; ++dd){
          robst[dd] = 0.;
        }
        // Reverse collision
        adjColOp(fdist[jj][ii],fdistOut,mac,obst[jj][ii],rfdistIn,
                 adj[jj][ii],rmac,robst);
        // Post processing
        for (size_t dd = 0; dd < dof; ++dd){
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
class NewAdjointThermalStreamingLoop2d : public AdjointLBLoop {

  using IsoLattice = typename CollisionOperator::IsoLatticeType;
  using ThermalLattice = typename CollisionOperator::ThermalLatticeType;
  using Real = typename CollisionOperator::RealType;
public:
  NewAdjointThermalStreamingLoop2d(Box2d _box, const Box2d& _boxGhosted)
    : boundingBox(_box)
  {
    xmin = _boxGhosted.xRange.getBeginId();
    xmax = _boxGhosted.xRange.getEndId();
    ymin = _boxGhosted.yRange.getBeginId();
    ymax = _boxGhosted.yRange.getEndId();
  }
  void execute(void* adj_v, void* adjPrev_v,void*,void*) override
  {
    auto adj = static_cast<Real***>(adj_v);
    auto adjPrev = static_cast<Real***>(adjPrev_v);
    PetscInt dd,nextx,nexty;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        // Iso streaming
        for (dd = 0; dd < IsoLattice::numDOF; ++dd){
          nextx = ii + IsoLattice::ex[dd];
          nexty = jj + IsoLattice::ey[dd];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              adj[jj][ii][dd] = adjPrev[nexty][nextx][dd];
            } else {
            adj[jj][ii][dd] = 0.;
          }
        }
        // Thermal streaming
        for (dd = 0; dd < ThermalLattice::numDOF; ++dd){
          nextx = ii + ThermalLattice::ex[dd];
          nexty = jj + ThermalLattice::ey[dd];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              adj[jj][ii][dd + IsoLattice::numDOF] =
                adjPrev[nexty][nextx][dd + IsoLattice::numDOF];
            } else {
            adj[jj][ii][dd + IsoLattice::numDOF] = 0.;
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
