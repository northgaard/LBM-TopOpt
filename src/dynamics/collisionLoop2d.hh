#ifndef COLLISIONLOOP2D
#define COLLISIONLOOP2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include <utility>

template <class CollisionOperator>
class CollisionLoop2d {

public:

  CollisionLoop2d(Box2d _box, CollisionOperator _op)
    : boundingBox(_box), colOp(_op){}
  void operator()(void* fdist_v, void* mac_v, void*){

    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;

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

template <class CollisionOperator>
class ObstCollisionLoop2d {

public:

  ObstCollisionLoop2d(Box2d _box, CollisionOperator _op)
    : boundingBox(_box), colOp(_op){}

  void operator()(void* fdist_v, void* mac_v, void* obst_v){

    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac;
    PetscScalar*** obst = (PetscScalar***) obst;

    for (auto jj : boundingBox.yRange){
      for (auto ii : boundingBox.xRange){
        colOp(fdist[jj][ii],fdist[jj][ii],mac[jj][ii],obst[jj][ii]);
      }
    }
  }

private:

  Box2d boundingBox;
  CollisionOperator colOp;

};

template <class Lattice, class CollisionOperator>
class CollisionAndSwapLoop2d {

public:

  CollisionAndSwapLoop2d(Box2d _box, CollisionOperator _op)
    : boundingBox(_box), colOp(_op){}

  void operator()(void* fdist_v, void* mac_v, void*){

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

  Box2d boundingBox;
  CollisionOperator colOp;

};

template <class Lattice, class CollisionOperator>
class ObstCollisionAndSwapLoop2d {

public:

  ObstCollisionAndSwapLoop2d(Box2d _box, CollisionOperator _op)
    : boundingBox(_box), colOp(_op){}

  void operator()(void* fdist_v, void* mac_v, void* obst_v){

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

  Box2d boundingBox;
  CollisionOperator colOp;

};

template <class Lattice>
class StreamBySwappingLoop2d {

public:

  StreamBySwappingLoop2d(Box2d _box) : boundingBox(_box){}

  void operator()(void* fdist_v,void*,void*){

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
class CollideAndStreamSingleLoop2d {

public:

  CollideAndStreamSingleLoop2d(Box2d _box, CollisionOperator _op)
    : boundingBox(_box), colOp(_op){}

  void operator()(void* fdist_v, void* mac_v, void*){

    PetscInt xmin = boundingBox.xRange.getBeginId();
    PetscInt xmax = boundingBox.xRange.getEndId();
    PetscInt ymin = boundingBox.yRange.getBeginId();
    PetscInt ymax = boundingBox.yRange.getEndId();
    PetscInt nextx,nexty,dd;

    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;
    PetscScalar ftemp;

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

  Box2d boundingBox;
  CollisionOperator colOp;

};

template <class Lattice, class CollisionOperator>
class ObstCollideAndStreamSingleLoop2d {

public:

  ObstCollideAndStreamSingleLoop2d(Box2d _box, CollisionOperator _op)
    : boundingBox(_box), colOp(_op){}

  void operator()(void* fdist_v, void* mac_v, void* obst_v){

    PetscInt xmin = boundingBox.xRange.getBeginId();
    PetscInt xmax = boundingBox.xRange.getEndId();
    PetscInt ymin = boundingBox.yRange.getBeginId();
    PetscInt ymax = boundingBox.yRange.getEndId();
    PetscInt nextx,nexty,dd;

    PetscScalar*** fdist = (PetscScalar***) fdist_v;
    PetscScalar*** mac = (PetscScalar***) mac_v;
    PetscScalar*** obst = (PetscScalar***) obst_v;
    PetscScalar ftemp;

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

  Box2d boundingBox;
  CollisionOperator colOp;

};

#endif
