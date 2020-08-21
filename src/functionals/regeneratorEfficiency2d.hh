#ifndef REGENERATOREFFICIENCY2D
#define REGENERATOREFFICIENCY2D

#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

template <class ThermalCollisionOperator>
class RegeneratorEfficiency2d : public LBMacroFunctional {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
  using Temperature = typename ThermalCollisionOperator::ThermalMacrosType;
public:
  void evaluate(PetscInt,const PetscScalar timeConst,
                void* mac_v, void* obst_v, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    PetscScalar pvn = 0.;
    PetscScalar pvd = 0.;
    PetscScalar vn = 0.;
    PetscScalar vd = 0.;

    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        pvn += timeConst*(1. - obst[jj][ii][0])*(mac[jj][ii][3] - Tinit);
        pvd += (1. - obst[jj][ii][0])*(Tblow - Tinit);
      }
    }
    MPI_Allreduce(&pvn,&vn,1,MPIU_SCALAR,MPI_SUM,communicator);
    MPI_Allreduce(&pvd,&vd,1,MPIU_SCALAR,MPI_SUM,communicator);
    *value = vn/vd;
  }
  void adjointCollideSource(PetscInt, const PetscScalar timeConst,
                            void* adj_v, void*, void* obst_v) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    // Compute denominator
    PetscScalar pvd = 0.;
    PetscScalar vd = 0.;

    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        pvd += timeConst*(1. - obst[jj][ii][0])*(Tblow - Tinit);
      }
    }
    MPI_Allreduce(&pvd,&vd,1,MPIU_SCALAR,MPI_SUM,communicator);
    // Add adjoint source
    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        for (size_t dd = 0; dd < ThermalLattice::numDOF; ++dd){
          adj[jj][ii][dd+IsoLattice::numDOF] += (1. - obst[jj][ii][0])/vd;
        }
      }
    }
  }
  void adjointSensitivitySource(PetscInt, const PetscScalar timeConst,
                                void* sens_v, void* fdist_v, void* obst_v) const override
  {
    auto sens = static_cast<PetscScalar**>(sens_v);
    auto fdist = static_cast<PetscScalar***>(fdist_v);
    auto obst = static_cast<PetscScalar***>(obst_v);
    PetscScalar pvn = 0.;
    PetscScalar pvd = 0.;
    PetscScalar vn = 0.;
    PetscScalar vd = 0.;
    PetscScalar T;
    // Compute nominator and denominator
    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        Temperature::compute(fdist[jj][ii] + IsoLattice::numDOF,T);
        pvn += timeConst*(1. - obst[jj][ii][0])*(T - Tinit);
        pvd += (1. - obst[jj][ii][0])*(Tblow - Tinit);
      }
    }
    MPI_Allreduce(&pvn,&vn,1,MPIU_SCALAR,MPI_SUM,communicator);
    MPI_Allreduce(&pvd,&vd,1,MPIU_SCALAR,MPI_SUM,communicator);
    // Add sensitivity source
    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        Temperature::compute(fdist[jj][ii] + IsoLattice::numDOF,T);
        sens[jj][ii] += ((Tblow - Tinit)*vn - (T - Tinit)*vd)/(vd*vd);
      }
    }
  }
  RegeneratorEfficiency2d(const Box2d globalDomain, const NewLBSolver2d& solver,
                          PetscScalar _Ti, PetscScalar _Tb)
    : Tinit(_Ti), Tblow(_Tb)
  {
    communicator = solver.communicator;
    boxIntersection(solver.getLocalBoundingBox(),globalDomain,domain);
  }
private:
  Box2d domain;
  PetscScalar Tinit,Tblow;
  MPI_Comm communicator;
};

#endif
