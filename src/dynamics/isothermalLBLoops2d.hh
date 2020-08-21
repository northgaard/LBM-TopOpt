#ifndef ISOTHERMALLBLOOPS2D
#define ISOTHERMALLBLOOPS2D

#include "petsc.h"
#include "LBLoop.hh"
#include "core/geometry2d.hh"
#include <utility>

template <class CollisionOperator>
class IsoCollisionLoop2d : public LBLoop {

public:
  IsoCollisionLoop2d(Box2d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto mac = static_cast<PetscScalar***>(mac_v);

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii]);
      }
    }
  }
private:
  CollisionOperator colOp;
  Box2d boundingBox;
};

template <class CollisionOperator>
class ObstIsoCollisionLoop2d : public LBLoop {

public:
  ObstIsoCollisionLoop2d(Box2d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto mac = static_cast<PetscScalar***>(mac_v);
    auto obst = static_cast<PetscScalar***>(obst_v);

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
      }
    }
  }
private:
  CollisionOperator colOp;
  Box2d boundingBox;
};

template <class Lattice>
class IsoStreamLoop2d : public LBLoop {

public:
  IsoStreamLoop2d(Box2d _box, const Box2d& _boxGlobal) : boundingBox(_box)
  {
    xmin = _boxGlobal.xRange.getBeginId();
    xmax = _boxGlobal.xRange.getEndId();
    ymin = _boxGlobal.yRange.getBeginId();
    ymax = _boxGlobal.yRange.getEndId();
  }
  void execute(void* fdist_g_v, void* fdist_l_v, void*) override
  {
    auto fdist_g = static_cast<PetscScalar***>(fdist_g_v);
    auto fdist_l = static_cast<const PetscScalar***>(fdist_l_v);
    PetscInt nextx,nexty,dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        for (dd = 1; dd < Lattice::numDOF; ++dd){
          nextx = ii - Lattice::ex[dd];
          nexty = jj - Lattice::ey[dd];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax)){
            fdist_g[jj][ii][dd] = fdist_l[nexty][nextx][dd];
          } else {
            // Default is mid-way bounceback
            fdist_g[jj][ii][dd] = fdist_l[jj][ii][Lattice::opp[dd]];
          }
        }
      }
    }
  }
private:
  Box2d boundingBox;
  PetscInt xmin,xmax,ymin,ymax;
};

template <class Lattice, class CollisionOperator>
class IsoCollisionAndSwapLoop2d : public LBLoop {

public:
  IsoCollisionAndSwapLoop2d(Box2d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;
    PetscInt dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii]);
        for (dd = 1; dd <= Lattice::half; ++dd){
          std::swap(fdist[jj][ii][dd],fdist[jj][ii][dd+Lattice::half]);
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box2d boundingBox;
};

template <class Lattice, class CollisionOperator>
class ObstIsoCollisionAndSwapLoop2d : public LBLoop {

public:
  ObstIsoCollisionAndSwapLoop2d(Box2d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;
    PetscScalar*** obst = (PetscScalar***) obst_v;
    PetscInt dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
        for (dd = 1; dd <= Lattice::half; ++dd){
          std::swap(fdist[jj][ii][dd],fdist[jj][ii][dd+Lattice::half]);
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box2d boundingBox;
};

template <class Lattice>
class IsoStreamBySwappingLoop2d : public LBLoop {

public:

  IsoStreamBySwappingLoop2d(Box2d _box) : boundingBox(_box){}
  void execute(void* fdist_v,void*,void*) override
  {
    PetscInt xmin = boundingBox.xRange.getBeginId();
    PetscInt xmax = boundingBox.xRange.getEndId();
    PetscInt ymin = boundingBox.yRange.getBeginId();
    PetscInt ymax = boundingBox.yRange.getEndId();
    PetscInt nextx,nexty,dd;
    PetscScalar*** fdist = (PetscScalar***) fdist_v;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        for (dd = 1; dd <= Lattice::half; ++dd){
          nextx = ii + Lattice::ex[dd];
          nexty = jj + Lattice::ey[dd];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax)){
            std::swap(fdist[jj][ii][dd+Lattice::half],fdist[nexty][nextx][dd]);
          }
        }
      }
    }
  }
private:
  Box2d boundingBox;
};

template <class Lattice, class CollisionOperator>
class IsoCollideAndStreamSingleLoop2d : public LBLoop {

public:

  IsoCollideAndStreamSingleLoop2d(Box2d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
  }
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;
    PetscScalar ftemp;
    PetscInt nextx,nexty,dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii]);
        for (dd = 1; dd <= Lattice::half; ++dd){
          nextx = ii + Lattice::ex[dd];
          nexty = jj + Lattice::ey[dd];
          ftemp = fdist[jj][ii][dd];
          fdist[jj][ii][dd] = fdist[jj][ii][dd+Lattice::half];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax)){
            fdist[jj][ii][dd+Lattice::half] = fdist[nexty][nextx][dd];
            fdist[nexty][nextx][dd] = ftemp;
          }
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box2d boundingBox;
  PetscInt xmin,xmax,ymin,ymax;
};

template <class Lattice, class CollisionOperator>
class ObstIsoCollideAndStreamSingleLoop2d : public LBLoop {

public:

  ObstIsoCollideAndStreamSingleLoop2d(Box2d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
  }
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;
    PetscScalar*** obst = (PetscScalar***) obst_v;
    PetscScalar ftemp;
    PetscInt nextx,nexty,dd;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
        for (dd = 1; dd <= Lattice::half; ++dd){
          nextx = ii + Lattice::ex[dd];
          nexty = jj + Lattice::ey[dd];
          ftemp = fdist[jj][ii][dd];
          fdist[jj][ii][dd] = fdist[jj][ii][dd+Lattice::half];
          if (PetscLikely(nextx >= xmin && nexty >= ymin &&
                          nextx <= xmax && nexty <= ymax)){
            fdist[jj][ii][dd+Lattice::half] = fdist[nexty][nextx][dd];
            fdist[nexty][nextx][dd] = ftemp;
          }
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box2d boundingBox;
  PetscInt xmin,xmax,ymin,ymax;
};

template <class Equilibrium>
class IsoEquilibriumInitializationLoop2d : public LBLoop {

public:
  IsoEquilibriumInitializationLoop2d(Box2d _box) : boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto mac = static_cast<PetscScalar***>(mac_v);
    PetscScalar uSqr;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        uSqr = mac[jj][ii][1]*mac[jj][ii][1] + mac[jj][ii][2]*mac[jj][ii][2];
        Equilibrium::setAllEquilibria(fdist[jj][ii],mac[jj][ii][0],mac[jj][ii][1],
                                      mac[jj][ii][2],uSqr);
      }
    }
  }
private:
  Box2d boundingBox;
};

#endif
