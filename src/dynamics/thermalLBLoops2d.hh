#ifndef THERMALLBLOOPS2D
#define THERMALLBLOOPS2D

#include "petsc.h"
#include "LBLoop.hh"
#include "core/geometry2d.hh"
#include <utility>
#include <type_traits>

template <class ThermalCollisionOperator>
class ObstThermalCollisionAndSwapLoop2d : public LBLoop {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  constexpr static bool thermalHasRestVelocity =
    (ThermalLattice::ex[0] == 0 && ThermalLattice::ey[0] == 0);
  constexpr static PetscInt thermalBegin =
    thermalHasRestVelocity ? 1 : 0;
  constexpr static PetscInt thermalEnd =
    (ThermalLattice::numDOF - 1) / 2;
public:
  ObstThermalCollisionAndSwapLoop2d(Box2d _box, ThermalCollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto mac = static_cast<PetscScalar***>(mac_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    PetscScalar* tdist;
    PetscInt dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
        for (dd = 1; dd <= IsoLattice::half; ++dd){
          std::swap(fdist[jj][ii][dd],fdist[jj][ii][dd+IsoLattice::half]);
        }
        tdist = fdist[jj][ii] + IsoLattice::numDOF;
        for (dd = thermalBegin; dd <= thermalEnd; ++dd){
          std::swap(tdist[dd],tdist[dd+ThermalLattice::half]);
        }
      }
    }
  }
private:
  ThermalCollisionOperator colOp;
  Box2d boundingBox;
};

template <class ThermalCollisionOperator>
class ThermalStreamBySwappingLoop2d : public LBLoop {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  constexpr static bool thermalHasRestVelocity =
    (ThermalLattice::ex[0] == 0 && ThermalLattice::ey[0] == 0);
  constexpr static PetscInt thermalBegin =
    thermalHasRestVelocity ? 1 : 0;
  constexpr static PetscInt thermalEnd =
    (ThermalLattice::numDOF - 1) / 2;
public:
  ThermalStreamBySwappingLoop2d(Box2d _box) : boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
  }
  void execute(void* fdist_v,void*,void*) override
  {
    PetscInt nextx,nexty,dd;
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscScalar *tdist, *tdistNext;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        for (dd = 1; dd <= IsoLattice::half; ++dd){
          nextx = ii + IsoLattice::ex[dd];
          nexty = jj + IsoLattice::ey[dd];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              std::swap(fdist[jj][ii][dd+IsoLattice::half],fdist[nexty][nextx][dd]);
            }
        }
        tdist = fdist[jj][ii] + IsoLattice::numDOF;
        for (dd = thermalBegin; dd <= thermalEnd; ++dd){
          nextx = ii + ThermalLattice::ex[dd];
          nexty = jj + ThermalLattice::ey[dd];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              tdistNext = fdist[nexty][nextx] + IsoLattice::numDOF;
              std::swap(tdist[dd+ThermalLattice::half],tdistNext[dd]);
            }
        }
      }
    }
  }
private:
  Box2d boundingBox;
  PetscInt xmin,xmax,ymin,ymax;
};

template <class ThermalCollisionOperator>
class ObstThermalCollideAndStreamSingleLoop2d : public LBLoop {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  constexpr static bool thermalHasRestVelocity =
    (ThermalLattice::ex[0] == 0 && ThermalLattice::ey[0] == 0);
  constexpr static PetscInt thermalBegin =
    thermalHasRestVelocity ? 1 : 0;
  constexpr static PetscInt thermalEnd =
    (ThermalLattice::numDOF - 1) / 2;
public:
  ObstThermalCollideAndStreamSingleLoop2d(Box2d _box, ThermalCollisionOperator _op)
    : colOp(_op), boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
  }
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto mac = static_cast<PetscScalar***>(mac_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    PetscScalar ftemp;
    PetscScalar *tdist, *tdistNext;
    PetscInt nextx,nexty,dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
        for (dd = 1; dd <= IsoLattice::half; ++dd){
          nextx = ii + IsoLattice::ex[dd];
          nexty = jj + IsoLattice::ey[dd];
          ftemp = fdist[jj][ii][dd];
          fdist[jj][ii][dd] = fdist[jj][ii][dd+IsoLattice::half];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              fdist[jj][ii][dd+IsoLattice::half] = fdist[nexty][nextx][dd];
              fdist[nexty][nextx][dd] = ftemp;
            }
        }
        tdist = fdist[jj][ii] + IsoLattice::numDOF;
        for (dd = thermalBegin; dd <= thermalEnd; ++dd){
          nextx = ii + ThermalLattice::ex[dd];
          nexty = jj + ThermalLattice::ey[dd];
          ftemp = tdist[dd];
          tdist[dd] = tdist[dd+ThermalLattice::half];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax))
            {
              tdistNext = fdist[nexty][nextx] + IsoLattice::numDOF;
              tdist[dd+ThermalLattice::half] = tdistNext[dd];
              tdistNext[dd] = ftemp;
            }
        }
      }
    }
  }
private:
  ThermalCollisionOperator colOp;
  Box2d boundingBox;
  PetscInt xmin,xmax,ymin,ymax;
};

template <class ThermalCollisionOperator>
class ThermalEquilibriumInitializationLoop2d : public LBLoop {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  using IsoEq = typename ThermalCollisionOperator::IsoEquilibriumType;
  using ThermalEq = typename ThermalCollisionOperator::ThermalEquilibriumType;
public:
  ThermalEquilibriumInitializationLoop2d(Box2d _box) : boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto mac = static_cast<PetscScalar***>(mac_v);
    // PetscInt dd;
    PetscScalar uSqr;
    PetscScalar *tdist;
    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        uSqr = mac[jj][ii][1]*mac[jj][ii][1] + mac[jj][ii][2]*mac[jj][ii][2];
        IsoEq::setAllEquilibria(fdist[jj][ii],mac[jj][ii][0],mac[jj][ii][1],mac[jj][ii][2],
                                uSqr);
        tdist = fdist[jj][ii] + IsoLattice::numDOF;
        ThermalEq::setAllEquilibria(tdist,mac[jj][ii][1],mac[jj][ii][2],mac[jj][ii][3]);
      }
    }
  }
private:
  Box2d boundingBox;
};

#endif
