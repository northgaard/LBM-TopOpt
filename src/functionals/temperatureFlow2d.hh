#ifndef TEMPERATUREFLOW2D
#define TEMPERATUREFLOW2D

#include "petsc.h"
#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

template <class ThermalCollisionOperator>
class ThermalOutflowEast2d : public LBMacroFunctional {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  using IsoMacros = typename ThermalCollisionOperator::IsoMacrosType;
  using Temperature = typename ThermalCollisionOperator::ThermalMacrosType;
public:
  void evaluate(PetscInt, const PetscScalar timeConst,
                void* mac_v, void*, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar pout = 0.;

    for (auto jj : outlet.yRange){
      for (auto ii : outlet.xRange){
        pout -= timeConst*mac[jj][ii][1]*mac[jj][ii][3];
      }
    }
    MPI_Allreduce(&pout,value,1,MPIU_SCALAR,MPI_SUM,communicator);
  }
  void adjointCollideSource(PetscInt, const PetscScalar timeConst,
                            void* adj_v, void* fdist_v, void*) const override
  {
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto adj = static_cast<PetscScalar***>(adj_v);
    PetscScalar* tdist;
    PetscScalar* tadj;
    PetscScalar ux,T;
    for (auto jj : outlet.yRange){
      for (auto ii : outlet.xRange){
        tdist = fdist[jj][ii] + IsoLattice::numDOF;
        tadj = adj[jj][ii] + IsoLattice::numDOF;
        IsoMacros::computeUx(fdist[jj][ii],ux);
        Temperature::compute(tdist,T);
        for (PetscInt dd = 0; dd < IsoLattice::numDOF; ++dd){
          adj[jj][ii][dd] -= timeConst*IsoLattice::ex[dd]*T;
        }
        for (PetscInt dd = 0; dd < ThermalLattice::numDOF; ++dd){
          tadj[dd] -= timeConst*ux;
        }
      }
    }
  }
  ThermalOutflowEast2d(const Box2d globalOutlet, const NewLBSolver2d& solver)
  {
    communicator = solver.communicator;
    boxIntersection(solver.getLocalBoundingBox(),globalOutlet,outlet);
  }
private:
  Box2d outlet;
  MPI_Comm communicator;
};

#endif
