#ifndef DOMAINTEMPERATURE2D
#define DOMAINTEMPERATURE2D

#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"

template <class ThermalCollisionOperator>
class DomainTemperature2d : public LBMacroFunctional {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
public:
  void evaluate(PetscInt,const PetscScalar timeConst,
                void* mac_v, void*, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar pv = 0.;

    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        pv += timeConst*invDomainArea*mac[jj][ii][3];
      }
    }
    MPI_Allreduce(&pv,value,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
  }
  void adjointCollideSource(PetscInt, const PetscScalar timeConst,
                            void* adj_v, void*,void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
        for (size_t dd = 0; dd < ThermalLattice::numDOF; ++dd){
          adj[jj][ii][dd+IsoLattice::numDOF] += timeConst*invDomainArea;
        }
      }
    }
  }
  DomainTemperature2d(const Box2d globalDomain, const Box2d localBoundingBox)
  {
    boxIntersection(localBoundingBox,globalDomain,domain);
    invDomainArea = (PetscScalar) (1./(globalDomain.getNx()*globalDomain.getNy()));
  }
private:
  Box2d domain;
  PetscScalar invDomainArea;
};

#endif
