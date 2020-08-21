#ifndef DOMAINTEMPERATURE2D
#define DOMAINTEMPERATURE2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "objectiveFunctions/objectiveFunction2d.hh"
#include "latticeBoltzmann/LBSolverBase2d.hh"
#include <memory>

template <class IsoLattice, class ThermalLattice>
class DomainTemperature2d : public ObjectiveFunction2d {

public:

  virtual void
  evaluate(PetscScalar*** mac, PetscScalar& value) const
  {
    value = 0.;
    PetscScalar local = 0.;

    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
	local -= invDomainArea * mac[jj][ii][3];
      }
    }

    MPI_Allreduce(&local,&value,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
    
  }

  virtual void
  adjointStateSource(PetscInt,PetscScalar*** adj,
		     PetscScalar*** fdist) const
  {

    PetscScalar* tadj;

    for (auto jj : domain.yRange){
      for (auto ii : domain.xRange){
	
	tadj = adj[jj][ii] + IsoLattice::numDOF;

	for (size_t dd = 0; dd < ThermalLattice::numDOF; ++dd){
	  tadj[dd] -= invDomainArea * invTimesteps;
	}
      }
    }
  }

  static std::unique_ptr<ObjectiveFunction2d>
  make (const Box2d fullDomain, const LBSolverBase2d& solver,
	PetscInt numTimesteps)
  {
    return std::unique_ptr<ObjectiveFunction2d>
      (new DomainTemperature2d<IsoLattice,ThermalLattice>
       (fullDomain,solver,numTimesteps));
  }

private:

  DomainTemperature2d(const Box2d fullDomain, const LBSolverBase2d& solver,
		     PetscInt numTimesteps)
  {
    boxIntersection(solver.getLocalBoundingBox(),fullDomain,domain);

    invDomainArea = (PetscScalar) (1./(fullDomain.getNx()*fullDomain.getNy()));
    invTimesteps = (PetscScalar) (1./numTimesteps);
  }

  Box2d domain;
  PetscScalar invDomainArea, invTimesteps; 

};

#endif
