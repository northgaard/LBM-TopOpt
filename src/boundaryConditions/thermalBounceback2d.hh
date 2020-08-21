#ifndef THERMALBOUNCEBACK2D
#define THERMALBOUNCEBACK2D

#include "petsc.h"
#include "boundaryConditions/LBBoundaryLoop.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "adjointBoundaryConditions/adjointThermalBounceback2d.hh"
#include "core/geometry2d.hh"
#include "core/macros.hh"
#include "latticeBoltzmann/lattices2d.hh"

/******
       Base templates
******/

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorth2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouth2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackWest2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackEast2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthWest2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthEast2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthWest2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthEast2d : public LBBoundaryLoop {};

/******
       D2Q4 implementation
******/

template <class IsoLattice>
class ThermalBouncebackNorth2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackNorth2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","North"),
                                         boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist;

    for (auto ii : boundingBox.xRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[0] = tdist[2];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackNorth2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackSouth2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackSouth2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","South"),
                                         boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist;

    for (auto ii : boundingBox.xRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[2] = tdist[0];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackSouth2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackWest2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackWest2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","West"),
                                        boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[3] = tdist[1];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackWest2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackEast2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackEast2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","East"),
                                        boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[1] = tdist[3];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackEast2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackNorthWest2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackNorthWest2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","NorthWest"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[0] = tdist[2];
    tdist[3] = tdist[1];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackNorthWest2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackNorthEast2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackNorthEast2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","NorthEast"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[0] = tdist[2];
    tdist[1] = tdist[3];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackNorthEast2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackSouthWest2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackSouthWest2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","SouthWest"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[2] = tdist[0];
    tdist[3] = tdist[1];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackSouthWest2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackSouthEast2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalBouncebackSouthEast2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","SouthEast"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[1] = tdist[3];
    tdist[2] = tdist[0];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackSouthEast2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

/******
       D2Q5 implementation
******/

template <class IsoLattice>
class ThermalBouncebackNorth2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackNorth2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","North"),
                                         boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist;

    for (auto ii : boundingBox.xRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[1] = tdist[3];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackNorth2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackSouth2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackSouth2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","South"),
                                         boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist;

    for (auto ii : boundingBox.xRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[3] = tdist[1];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackSouth2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackWest2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackWest2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","West"),
                                        boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[4] = tdist[2];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackWest2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackEast2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackEast2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","East"),
                                        boundingBox(_box) {}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tdist[2] = tdist[4];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackEast2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackNorthWest2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackNorthWest2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","NorthWest"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[1] = tdist[3];
    tdist[4] = tdist[2];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackNorthWest2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackNorthEast2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackNorthEast2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","NorthEast"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[1] = tdist[3];
    tdist[2] = tdist[4];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackNorthEast2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackSouthWest2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackSouthWest2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","SouthWest"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[3] = tdist[1];
    tdist[4] = tdist[2];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackSouthWest2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalBouncebackSouthEast2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalBouncebackSouthEast2d(Box2d _box) : LBBoundaryLoop("ThermalBounceback2d","SouthEast"),
                                             boundingBox(_box){}
  void execute(PetscInt, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscInt jj = boundingBox.yRange.getBeginId();
    PetscScalar* tdist = fdist[jj][ii] + IsoLattice::numDOF;
    tdist[2] = tdist[4];
    tdist[3] = tdist[1];
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalBouncebackSouthEast2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
};

/******
       Descriptor
******/

template <class IsoLattice, class ThermalLattice>
class ThermalBounceback2d : public NewBoundaryDescriptor2d {};

template <class IsoLattice>
class ThermalBounceback2d<IsoLattice,D2Q4> : public NewBoundaryDescriptor2d {

public:
  ThermalBounceback2d() : NewBoundaryDescriptor2d("ThermalBouncebac2d"){}
  PetscErrorCode getNorthBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackNorth2d<IsoLattice,D2Q4>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackSouth2d<IsoLattice,D2Q4>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackWest2d<IsoLattice,D2Q4>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackEast2d<IsoLattice,D2Q4>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthWestBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackNorthWest2d<IsoLattice,D2Q4>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthEastBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackNorthEast2d<IsoLattice,D2Q4>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthWestBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackSouthWest2d<IsoLattice,D2Q4>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthEastBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackSouthEast2d<IsoLattice,D2Q4>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
};

template <class IsoLattice>
class ThermalBounceback2d<IsoLattice,D2Q5> : public NewBoundaryDescriptor2d {

public:
  ThermalBounceback2d() : NewBoundaryDescriptor2d("ThermalBouncebac2d"){}
  PetscErrorCode getNorthBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackNorth2d<IsoLattice,D2Q5>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthBoundary(const Box2d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackSouth2d<IsoLattice,D2Q5>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackWest2d<IsoLattice,D2Q5>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackEast2d<IsoLattice,D2Q5>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthWestBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackNorthWest2d<IsoLattice,D2Q5>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthEastBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackNorthEast2d<IsoLattice,D2Q5>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthWestBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackSouthWest2d<IsoLattice,D2Q5>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthEastBoundary(const Box2d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalBouncebackSouthEast2d<IsoLattice,D2Q5>
      (boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
};

#endif
