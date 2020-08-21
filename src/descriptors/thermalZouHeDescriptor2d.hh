#ifndef THERMALZOUHEDESCRIPTOR2D
#define THERMALZOUHEDESCRIPTOR2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "boundaryConditions/thermalZouHe2d.hh"
#include "adjointBoundaryConditions/adjointThermalZouHe2d.hh"
#include <memory>

template <class IsoLattice, class ThermalLattice, class Input>
class ThermalZouHe2d : public BoundaryDescriptor2d {

public:

  static std::unique_ptr<BoundaryDescriptor2d> make(Input _tf)
  {
    return std::unique_ptr<BoundaryDescriptor2d>
      (new ThermalZouHe2d<IsoLattice,ThermalLattice,Input>(_tf));
  }

  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
				  BoundaryFunction& func) const override
  {

    ThermalZouHeWest2d<IsoLattice,ThermalLattice,Input>
      tzhWest(boundaryLocation,tempFunctional);
    func = [tzhWest](PetscInt t, PetscScalar*** fdist)
      { tzhWest(t,fdist); };
    return 0;

  }

  PetscErrorCode getAdjointWestBoundary(const Box2d boundaryLocation,
					AdjointBoundaryFunction& func) const override
  {

    AdjointThermalZouHeWest2d<IsoLattice,ThermalLattice>
      atzhWest(boundaryLocation);
    func = [atzhWest](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atzhWest(t,adj,fdist); };
    return 0;

  }

private:

  ThermalZouHe2d(Input _tf)
    : BoundaryDescriptor2d("ThermalZouHe2d"), tempFunctional(_tf) {}
  Input tempFunctional;

};

template <class IsoLattice, class ThermalLattice>
class ThermalZouHeNeumann2d : public BoundaryDescriptor2d {

public:

  static std::unique_ptr<BoundaryDescriptor2d> make()
  {
    return std::unique_ptr<BoundaryDescriptor2d>
      (new ThermalZouHeNeumann2d<IsoLattice,ThermalLattice>());
  }

  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
				 BoundaryFunction& func) const override
  {
    
    ThermalZouHeNeumannEast2d<IsoLattice,ThermalLattice>
      tzhnEast(boundaryLocation);
    func = [tzhnEast](PetscInt t, PetscScalar*** fdist)
      { tzhnEast(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointEastBoundary(const Box2d boundaryLocation,
					AdjointBoundaryFunction& func) const override
  {
    AdjointThermalZouHeNeumannEast2d<IsoLattice,ThermalLattice>
      atzhnEast(boundaryLocation);
    func = [atzhnEast](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atzhnEast(t,adj,fdist); };
    return 0;
    
  }

  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
				 BoundaryFunction& func) const override
  {

    ThermalZouHeNeumannWest2d<IsoLattice,ThermalLattice>
      tzhnWest(boundaryLocation);
    func = [tzhnWest](PetscInt t, PetscScalar*** fdist)
      { tzhnWest(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointWestBoundary(const Box2d boundaryLocation,
					AdjointBoundaryFunction& func) const override
  {
    
    AdjointThermalZouHeNeumannWest2d<IsoLattice,ThermalLattice>
      atzhnWest(boundaryLocation);
    func = [atzhnWest](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atzhnWest(t,adj,fdist); };
    return 0;
    
  }

private:

  ThermalZouHeNeumann2d() : BoundaryDescriptor2d("ThermalZouHeNeumann2d") {}
  
};

#endif
