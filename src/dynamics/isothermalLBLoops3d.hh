#ifndef ISOTHERMALLBLOOPS3D
#define ISOTHERMALLBLOOPS3D

#include "petsc.h"
#include "LBLoop.hh"
#include "core/geometry3d.hh"
#include "traits/collisionTraits3d.hh"
#include <utility>

template <class CollisionOperator>
class IsoCollisionLoop3d : public LBLoop {

public:
  IsoCollisionLoop3d(Box3d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          colOp(fdist[kk][jj][ii],fdist[kk][jj][ii],mac[kk][jj][ii]);
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box3d boundingBox;
};

template <class CollisionOperator>
class ObstIsoCollisionLoop3d : public LBLoop {

public:
  ObstIsoCollisionLoop3d(Box3d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);
    auto obst = static_cast<PetscScalar****>(obst_v);

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          colOp(fdist[kk][jj][ii],fdist[kk][jj][ii],
                mac[kk][jj][ii],obst[kk][jj][ii]);
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box3d boundingBox;
};

template <class CollisionOperator>
class IsoCollisionAndSwapLoop3d : public LBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
public:
  IsoCollisionAndSwapLoop3d(Box3d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);
    PetscInt dd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          colOp(fdist[kk][jj][ii],fdist[kk][jj][ii],mac[kk][jj][ii]);
          for (dd = 1; dd <= Lattice::half; ++dd){
            std::swap(fdist[kk][jj][ii][dd],fdist[kk][jj][ii][dd+Lattice::half]);
          }
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box3d boundingBox;
};

template <class CollisionOperator>
class ObstIsoCollisionAndSwapLoop3d : public LBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
public:
  ObstIsoCollisionAndSwapLoop3d(Box3d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);
    auto obst = static_cast<PetscScalar****>(obst_v);
    PetscInt dd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          colOp(fdist[kk][jj][ii],fdist[kk][jj][ii],
                mac[kk][jj][ii],obst[kk][jj][ii]);
          for (dd = 1; dd <= Lattice::half; ++dd){
            std::swap(fdist[kk][jj][ii][dd],fdist[kk][jj][ii][dd+Lattice::half]);
          }
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box3d boundingBox;
};

template <class CollisionOperator>
class IsoStreamBySwappingLoop3d : public LBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
public:
  IsoStreamBySwappingLoop3d(Box3d _box) : boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
    zmin = boundingBox.zRange.getBeginId();
    zmax = boundingBox.zRange.getEndId();
  }
  void execute(void* fdist_v,void*,void*) override
  {
    PetscInt nextx,nexty,nextz,dd;
    auto fdist = static_cast<PetscScalar****>(fdist_v);

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          for (dd = 1; dd <= Lattice::half; ++dd){
            nextx = ii + Lattice::ex[dd];
            nexty = jj + Lattice::ey[dd];
            nextz = kk + Lattice::ez[dd];
            if (PetscLikely(nextx >= xmin && nexty >= ymin && nextz >= zmin &&
                            nextx <= xmax && nexty <= ymax && nextz <= zmax))
              {
                std::swap(fdist[kk][jj][ii][dd+Lattice::half],
                          fdist[nextz][nexty][nextx][dd]);
              }
          }
        }
      }
    }
  }
private:
  Box3d boundingBox;
  PetscInt xmin,xmax,ymin,ymax,zmin,zmax;
};

template <class CollisionOperator>
class IsoCollideAndStreamSingleLoop3d : public LBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
public:
  IsoCollideAndStreamSingleLoop3d(Box3d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
    zmin = boundingBox.zRange.getBeginId();
    zmax = boundingBox.zRange.getEndId();
  }
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);
    PetscScalar ftemp;
    PetscInt nextx,nexty,nextz,dd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          colOp(fdist[kk][jj][ii],fdist[kk][jj][ii],mac[kk][jj][ii]);
          for (dd = 1; dd <= Lattice::half; ++dd){
            nextx = ii + Lattice::ex[dd];
            nexty = jj + Lattice::ey[dd];
            nextz = kk + Lattice::ez[dd];
            ftemp = fdist[kk][jj][ii][dd];
            fdist[kk][jj][ii][dd] = fdist[kk][jj][ii][dd+Lattice::half];
            if (PetscLikely(nextx >= xmin && nexty >= ymin && nextz >= zmin &&
                            nextx <= xmax && nexty <= ymax && nextz <= zmax))
              {
                fdist[kk][jj][ii][dd+Lattice::half] =
                  fdist[nextz][nexty][nextx][dd];
                fdist[nextz][nexty][nextx][dd] = ftemp;
              }
          }
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box3d boundingBox;
  PetscInt xmin,xmax,ymin,ymax,zmin,zmax;
};

template <class CollisionOperator>
class ObstIsoCollideAndStreamSingleLoop3d : public LBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
public:
  ObstIsoCollideAndStreamSingleLoop3d(Box3d _box, CollisionOperator _op)
    : colOp(_op), boundingBox(_box)
  {
    xmin = boundingBox.xRange.getBeginId();
    xmax = boundingBox.xRange.getEndId();
    ymin = boundingBox.yRange.getBeginId();
    ymax = boundingBox.yRange.getEndId();
    zmin = boundingBox.zRange.getBeginId();
    zmax = boundingBox.zRange.getEndId();
  }
  void execute(void* fdist_v, void* mac_v, void* obst_v) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);
    auto obst = static_cast<PetscScalar****>(obst_v);
    PetscScalar ftemp;
    PetscInt nextx,nexty,nextz,dd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          colOp(fdist[kk][jj][ii],fdist[kk][jj][ii],
                mac[kk][jj][ii],obst[kk][jj][ii]);
          for (dd = 1; dd <= Lattice::half; ++dd){
            nextx = ii + Lattice::ex[dd];
            nexty = jj + Lattice::ey[dd];
            nextz = kk + Lattice::ez[dd];
            ftemp = fdist[kk][jj][ii][dd];
            fdist[kk][jj][ii][dd] = fdist[kk][jj][ii][dd+Lattice::half];
            if (PetscLikely(nextx >= xmin && nexty >= ymin && nextz >= zmin &&
                            nextx <= xmax && nexty <= ymax && nextz <= zmax))
              {
                fdist[kk][jj][ii][dd+Lattice::half] =
                  fdist[nextz][nexty][nextx][dd];
                fdist[nextz][nexty][nextx][dd] = ftemp;
              }
          }
        }
      }
    }
  }
private:
  CollisionOperator colOp;
  Box3d boundingBox;
  PetscInt xmin,xmax,ymin,ymax,zmin,zmax;
};

template <class CollisionOperator>
class IsoEquilibriumInitializationLoop3d : public LBLoop {

  using Lattice = typename get_lattice<CollisionOperator>::type;
  using Equilibrium = typename CollisionOperator::Equilibrium;
public:
  IsoEquilibriumInitializationLoop3d(Box3d _box) : boundingBox(_box){}
  void execute(void* fdist_v, void* mac_v, void*) override
  {
    auto fdist = static_cast<PetscScalar****>(fdist_v);
    auto mac = static_cast<PetscScalar****>(mac_v);
    PetscScalar uSqr;
    PetscInt dd;

    for (auto kk : boundingBox.zRange){
      for (auto jj : boundingBox.yRange){
        for (auto ii : boundingBox.xRange){
          uSqr = mac[kk][jj][ii][1]*mac[kk][jj][ii][1]+
            mac[kk][jj][ii][2]*mac[kk][jj][ii][2]+
            mac[kk][jj][ii][3]*mac[kk][jj][ii][3];
          for (dd = 0; dd < Lattice::numDOF; ++dd){
            fdist[kk][jj][ii][dd] = Equilibrium::eq(dd,mac[kk][jj][ii][0],mac[kk][jj][ii][1],
                                                    mac[kk][jj][ii][2],mac[kk][jj][ii][3],uSqr);
          }
        }
      }
    }
  }
private:
  Box3d boundingBox;
};

#endif
