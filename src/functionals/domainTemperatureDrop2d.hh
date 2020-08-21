#ifndef DOMAINTEMPERATUREDROP2D
#define DOMAINTEMPERATUREDROP2D

#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"

template <class ThermalCollisionOperator>
class DomainTemperatureDrop2d : public LBMacroFunctional {

  using IsoLattice = typename ThermalCollisionOperator::IsoLatticeType;
  using ThermalLattice = typename ThermalCollisionOperator::ThermalLatticeType;
public:
  void evaluate(PetscInt, const PetscScalar timeConst,
                void* mac_v, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar pv = 0.;

    for (auto jj : domainIn.yRange){
      for (auto ii : domainIn.xRange){
        pv += timeConst*invDomainInArea*mac[jj][ii][3];
      }
    }
    for (auto jj : domainOut.yRange){
      for (auto ii : domainOut.xRange){
        pv -= timeConst*invDomainOutArea*mac[jj][ii][3];
      }
    }
    MPI_Allreduce(&pv,value,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
  }
  void adjointCollideSource(PetscInt, const PetscScalar timeConst,
                            void* adj_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    for (auto jj : domainIn.yRange){
      for (auto ii : domainIn.xRange){
        for (size_t dd = 0; dd < ThermalLattice::numDOF; ++dd){
          adj[jj][ii][dd+IsoLattice::numDOF] += timeConst*invDomainInArea;
        }
      }
    }
    for (auto jj : domainOut.yRange){
      for (auto ii : domainOut.xRange){
        for (size_t dd = 0; dd < ThermalLattice::numDOF; ++dd){
          adj[jj][ii][dd+IsoLattice::numDOF] -= timeConst*invDomainOutArea;
        }
      }
    }
  }
  DomainTemperatureDrop2d(const Box2d globalDomainIn, const Box2d globalDomainOut,
                          const Box2d localBoundingBox)
  {
    boxIntersection(localBoundingBox,globalDomainIn,domainIn);
    boxIntersection(localBoundingBox,globalDomainOut,domainOut);
    invDomainInArea = (PetscScalar) (1./(globalDomainIn.getNx()*globalDomainIn.getNy()));
    invDomainOutArea = (PetscScalar) (1./(globalDomainOut.getNx()*globalDomainOut.getNy()));
  }
private:
  Box2d domainIn;
  Box2d domainOut;
  PetscScalar invDomainInArea;
  PetscScalar invDomainOutArea;
};

#endif
