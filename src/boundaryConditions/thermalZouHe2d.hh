#ifndef THERMALZOUHE2D
#define THERMALZOUHE2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "boundaryConditions/LBBoundaryLoop.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "adjointBoundaryConditions/adjointThermalZouHe2d.hh"
#include "dynamics/thermalMacros2d.hh"
#include "latticeBoltzmann/lattices2d.hh"
#include "core/macros.hh"

/************ Zou He standard boundaries ************/

/******
       Base templates
******/

template <class IsoLattice, class ThermalLattice, class Input>
class ThermalZouHeWest2d : public LBBoundaryLoop {};

/******
       D2Q4 implementation
******/

template <class IsoLattice, class Input>
class ThermalZouHeWest2d<IsoLattice,D2Q4,Input> : public LBBoundaryLoop {

public:
  ThermalZouHeWest2d(Box2d _box, Input _temp)
    : LBBoundaryLoop("ThermalZouHe2d","West"), boundingBox(_box), tempFunc(_temp) {}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar T;
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tempFunc(timestep,ii,jj,T);
      tdist[3] = T - tdist[0] - tdist[1] - tdist[2];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalZouHeWest2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input tempFunc;
};

/******
       D2Q5 implementation
******/

template <class IsoLattice, class Input>
class ThermalZouHeWest2d<IsoLattice,D2Q5,Input> : public LBBoundaryLoop {

public:
  ThermalZouHeWest2d(Box2d _box, Input _temp)
    : LBBoundaryLoop("ThermalZouHe2d","West"), boundingBox(_box), tempFunc(_temp) {}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar T;
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      tempFunc(timestep,ii,jj,T);
      tdist[4] = T - tdist[0] - tdist[1] - tdist[2] - tdist[3];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalZouHeWest2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  Box2d boundingBox;
  Input tempFunc;
};

/******
       Boundary descriptor
******/

template <class IsoLattice, class ThermalLattice, class Input>
class ThermalZouHe2d : public NewBoundaryDescriptor2d {};

template <class IsoLattice, class Input>
class ThermalZouHe2d<IsoLattice,D2Q4,Input> : public NewBoundaryDescriptor2d {

public:
  ThermalZouHe2d(Input _tf) : NewBoundaryDescriptor2d("ThermalZouHe2d"), tempFunc(_tf){}
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalZouHeWest2d<IsoLattice,D2Q4,Input>
      (boundaryLocation,tempFunc);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
private:
  Input tempFunc;
};

template <class IsoLattice, class Input>
class ThermalZouHe2d<IsoLattice,D2Q5,Input> : public NewBoundaryDescriptor2d {

public:
  ThermalZouHe2d(Input _tf) : NewBoundaryDescriptor2d("ThermalZouHe2d"), tempFunc(_tf){}
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalZouHeWest2d<IsoLattice,D2Q5,Input>
      (boundaryLocation,tempFunc);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
private:
  Input tempFunc;
};

/************ Zou He Neumann boundaries ************/

/******
       Base templates
******/

template <class IsoLattice, class ThermalLattice>
class ThermalZouHeNeumannEast2d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalZouHeNeumannWest2d : public LBBoundaryLoop {};

/******
       D2Q4 implementation
******/

template <class IsoLattice>
class ThermalZouHeNeumannEast2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalZouHeNeumannEast2d(Box2d _box) : LBBoundaryLoop("ThermalZouHeNeumann2d","East"),
                                          boundingBox(_box) {}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar T,Tp1 = 0.,Tp2 = 0.;
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii-1] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp1);
      tdist = fdist[jj][ii-2] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp2);
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      T = (4.*Tp1 - Tp2)/3.;
      tdist[1] = T - tdist[0] - tdist[2] - tdist[3];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalZouHeNeumannEast2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  using Macros = StandardThermalMacros2d<D2Q4>;
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalZouHeNeumannWest2d<IsoLattice,D2Q4> : public LBBoundaryLoop {

public:
  ThermalZouHeNeumannWest2d(Box2d _box) : LBBoundaryLoop("ThermalZouHeNeumann2d","West"),
                                          boundingBox(_box) {}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar T,Tp1 = 0.,Tp2 = 0.;
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii+1] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp1);
      tdist = fdist[jj][ii+2] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp2);
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      T = (4.*Tp1 - Tp2)/3.;
      tdist[3] = T - tdist[0] - tdist[1] - tdist[2];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalZouHeNeumannWest2d<IsoLattice,D2Q4>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  using Macros = StandardThermalMacros2d<D2Q4>;
  Box2d boundingBox;
};

/******
       D2Q5 implementation
******/

template <class IsoLattice>
class ThermalZouHeNeumannEast2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalZouHeNeumannEast2d(Box2d _box) : LBBoundaryLoop("ThermalZouHeNeumann2d","East"),
                                          boundingBox(_box) {}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar T,Tp1 = 0.,Tp2 = 0.;
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii-1] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp1);
      tdist = fdist[jj][ii-2] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp2);
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      T = (4.*Tp1 - Tp2)/3.;
      tdist[2] = T - tdist[0] - tdist[1] - tdist[3] - tdist[4];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalZouHeNeumannEast2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  using Macros = StandardThermalMacros2d<D2Q5>;
  Box2d boundingBox;
};

template <class IsoLattice>
class ThermalZouHeNeumannWest2d<IsoLattice,D2Q5> : public LBBoundaryLoop {

public:
  ThermalZouHeNeumannWest2d(Box2d _box) : LBBoundaryLoop("ThermalZouHeNeumann2d","West"),
                                          boundingBox(_box) {}
  void execute(PetscInt timestep, void* fdist_v) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    PetscInt ii = boundingBox.xRange.getBeginId();
    PetscScalar T,Tp1 = 0.,Tp2 = 0.;
    PetscScalar* tdist;

    for (auto jj : boundingBox.yRange){
      tdist = fdist[jj][ii+1] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp1);
      tdist = fdist[jj][ii+2] + IsoLattice::numDOF;
      Macros::compute(tdist,Tp2);
      tdist = fdist[jj][ii] + IsoLattice::numDOF;
      T = (4.*Tp1 - Tp2)/3.;
      tdist[4] = T - tdist[0] - tdist[1] - tdist[2] - tdist[3];
    }
  }
  PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const override
  {
    PetscFunctionBeginUser;
    *adjLoop = new (std::nothrow) AdjointThermalZouHeNeumannWest2d<IsoLattice,D2Q5>
      (boundingBox);
    CHKNEWPTR(*adjLoop);
    PetscFunctionReturn(0);
  }
private:
  using Macros = StandardThermalMacros2d<D2Q5>;
  Box2d boundingBox;
};

/******
       Boundary descriptor
******/

template <class IsoLattice, class ThermalLattice>
class ThermalZouHeNeumann2d : public NewBoundaryDescriptor2d {};

template <class IsoLattice>
class ThermalZouHeNeumann2d<IsoLattice,D2Q4> : public NewBoundaryDescriptor2d {

public:
  ThermalZouHeNeumann2d() : NewBoundaryDescriptor2d("ThermalZouHeNeumann2d"){}
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalZouHeNeumannWest2d<IsoLattice,D2Q4>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalZouHeNeumannEast2d<IsoLattice,D2Q4>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
};

template <class IsoLattice>
class ThermalZouHeNeumann2d<IsoLattice,D2Q5> : public NewBoundaryDescriptor2d {

public:
  ThermalZouHeNeumann2d() : NewBoundaryDescriptor2d("ThermalZouHeNeumann2d"){}
  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalZouHeNeumannWest2d<IsoLattice,D2Q5>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) ThermalZouHeNeumannEast2d<IsoLattice,D2Q5>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
};

#endif
