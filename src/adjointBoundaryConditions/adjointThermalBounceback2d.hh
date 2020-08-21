#ifndef ADJOINTTHERMALBOUNCEBACK2D
#define ADJOINTTHERMALBOUNCEBACK2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"
#include "latticeBoltzmann/lattices2d.hh"

/******
       Base templates
******/

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackNorth2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackSouth2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackWest2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackEast2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackNorthWest2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackNorthEast2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackSouthWest2d : public AdjointLBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class AdjointThermalBouncebackSouthEast2d : public AdjointLBBoundaryLoop {};

/******
       D2Q4 implementation
******/

template <class IsoLattice>
class AdjointThermalBouncebackNorth2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackNorth2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj;

    for (auto ii : boundingBox.xRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[2] += tadj[0];
      tadj[0] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackSouth2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackSouth2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj;

    for (auto ii : boundingBox.xRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[0] += tadj[2];
      tadj[2] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackWest2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackWest2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tadj;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[1] += tadj[3];
      tadj[3] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackEast2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackEast2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tadj;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[3] += tadj[1];
      tadj[1] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackNorthWest2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackNorthWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[2] += tadj[0];
    tadj[1] += tadj[3];
    tadj[0] = 0.;
    tadj[3] = 0.;
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackNorthEast2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackNorthEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[2] += tadj[0];
    tadj[3] += tadj[1];
    tadj[0] = 0.;
    tadj[1] = 0.;
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackSouthWest2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackSouthWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[0] += tadj[2];
    tadj[1] += tadj[3];
    tadj[2] = 0.;
    tadj[3] = 0.;
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackSouthEast2d<IsoLattice,D2Q4> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackSouthEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[3] += tadj[1];
    tadj[0] += tadj[2];
    tadj[1] = 0.;
    tadj[2] = 0.;
  }
private:
  Box2d boundingBox;
};

/******
       D2Q5 implementation
******/

template <class IsoLattice>
class AdjointThermalBouncebackNorth2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackNorth2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj;

    for (auto ii : boundingBox.xRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[3] += tadj[1];
      tadj[1] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackSouth2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackSouth2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj;

    for (auto ii : boundingBox.xRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[1] += tadj[3];
      tadj[3] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackWest2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackWest2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tadj;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[2] += tadj[4];
      tadj[4] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackEast2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackEast2d(Box2d _box) : boundingBox(_box) {}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tadj;

    for (auto jj : boundingBox.yRange){
      tadj = adj[jj][ii] + IsoLattice::numDOF;
      tadj[4] += tadj[2];
      tadj[2] = 0.;
    }
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackNorthWest2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackNorthWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[3] += tadj[1];
    tadj[2] += tadj[4];
    tadj[1] = 0.;
    tadj[4] = 0.;
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackNorthEast2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackNorthEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[3] += tadj[1];
    tadj[4] += tadj[2];
    tadj[1] = 0.;
    tadj[2] = 0.;
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackSouthWest2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackSouthWest2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[1] += tadj[3];
    tadj[2] += tadj[4];
    tadj[3] = 0.;
    tadj[4] = 0.;
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class AdjointThermalBouncebackSouthEast2d<IsoLattice,D2Q5> : public AdjointLBBoundaryLoop {

public:
  AdjointThermalBouncebackSouthEast2d(Box2d _box) : boundingBox(_box){}
  void execute(PetscInt, void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tadj = adj[jj][ii] + IsoLattice::numDOF;
    tadj[4] += tadj[2];
    tadj[1] += tadj[3];
    tadj[2] = 0.;
    tadj[3] = 0.;
  }
private:
  Box2d boundingBox;
};

#endif
