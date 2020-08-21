#ifndef THERMALBOUNCEBACKDESCRIPTOR2D
#define THERMALBOUNCEBACKDESCRIPTOR2D

#include "petsc.h"
#include "core/geometry2d.hh"
#include "boundaryConditions/boundaryDescriptor2d.hh"
#include "boundaryConditions/thermalBounceback2d.hh"
#include "adjointBoundaryConditions/adjointThermalBounceback2d.hh"
#include <memory>

template <class IsoLattice, class ThermalLattice>
class ThermalBounceback2d : public BoundaryDescriptor2d {

public:

  static std::unique_ptr<BoundaryDescriptor2d> make()
  {
    return std::unique_ptr<BoundaryDescriptor2d>
      (new ThermalBounceback2d<IsoLattice,ThermalLattice>());
  }

  PetscErrorCode getNorthBoundary(const Box2d boundaryLocation,
				  BoundaryFunction& func) const override
  {

    ThermalBouncebackNorth2d<IsoLattice,ThermalLattice>
      tbbNorth(boundaryLocation);
    func = [tbbNorth](PetscInt t, PetscScalar*** fdist)
      { tbbNorth(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointNorthBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackNorth2d<IsoLattice,ThermalLattice>
      atbbNorth(boundaryLocation);
    func = [atbbNorth](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbNorth(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getSouthBoundary(const Box2d boundaryLocation,
				  BoundaryFunction& func) const override
  {

    ThermalBouncebackSouth2d<IsoLattice,ThermalLattice>
      tbbSouth(boundaryLocation);
    func = [tbbSouth](PetscInt t, PetscScalar*** fdist)
      { tbbSouth(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointSouthBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackSouth2d<IsoLattice,ThermalLattice>
      atbbSouth(boundaryLocation);
    func = [atbbSouth](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbSouth(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getWestBoundary(const Box2d boundaryLocation,
				  BoundaryFunction& func) const override
  {

    ThermalBouncebackWest2d<IsoLattice,ThermalLattice>
      tbbWest(boundaryLocation);
    func = [tbbWest](PetscInt t, PetscScalar*** fdist)
      { tbbWest(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointWestBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackWest2d<IsoLattice,ThermalLattice>
      atbbWest(boundaryLocation);
    func = [atbbWest](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbWest(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getEastBoundary(const Box2d boundaryLocation,
				  BoundaryFunction& func) const override
  {

    ThermalBouncebackEast2d<IsoLattice,ThermalLattice>
      tbbEast(boundaryLocation);
    func = [tbbEast](PetscInt t, PetscScalar*** fdist)
      { tbbEast(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointEastBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackEast2d<IsoLattice,ThermalLattice>
      atbbEast(boundaryLocation);
    func = [atbbEast](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbEast(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getNorthWestBoundary(const Box2d boundaryLocation,
				      BoundaryFunction& func) const override
  {

    ThermalBouncebackNorthWest2d<IsoLattice,ThermalLattice>
      tbbNorthWest(boundaryLocation);
    func = [tbbNorthWest](PetscInt t, PetscScalar*** fdist)
      { tbbNorthWest(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointNorthWestBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackNorthWest2d<IsoLattice,ThermalLattice>
      atbbNorthWest(boundaryLocation);
    func = [atbbNorthWest](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbNorthWest(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getNorthEastBoundary(const Box2d boundaryLocation,
				      BoundaryFunction& func) const override
  {

    ThermalBouncebackNorthEast2d<IsoLattice,ThermalLattice>
      tbbNorthEast(boundaryLocation);
    func = [tbbNorthEast](PetscInt t, PetscScalar*** fdist)
      { tbbNorthEast(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointNorthEastBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackNorthEast2d<IsoLattice,ThermalLattice>
      atbbNorthEast(boundaryLocation);
    func = [atbbNorthEast](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbNorthEast(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getSouthWestBoundary(const Box2d boundaryLocation,
				      BoundaryFunction& func) const override
  {

    ThermalBouncebackSouthWest2d<IsoLattice,ThermalLattice>
      tbbSouthWest(boundaryLocation);
    func = [tbbSouthWest](PetscInt t, PetscScalar*** fdist)
      { tbbSouthWest(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointSouthWestBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackSouthWest2d<IsoLattice,ThermalLattice>
      atbbSouthWest(boundaryLocation);
    func = [atbbSouthWest](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbSouthWest(t,adj,fdist); };
    return 0;

  }

  PetscErrorCode getSouthEastBoundary(const Box2d boundaryLocation,
				      BoundaryFunction& func) const override
  {

    ThermalBouncebackSouthEast2d<IsoLattice,ThermalLattice>
      tbbSouthEast(boundaryLocation);
    func = [tbbSouthEast](PetscInt t, PetscScalar*** fdist)
      { tbbSouthEast(t,fdist); };
    return 0;
    
  }

  PetscErrorCode getAdjointSouthEastBoundary(const Box2d boundaryLocation,
					 AdjointBoundaryFunction& func) const override
  {

    AdjointThermalBouncebackSouthEast2d<IsoLattice,ThermalLattice>
      atbbSouthEast(boundaryLocation);
    func = [atbbSouthEast](PetscInt t, PetscScalar*** adj, PetscScalar*** fdist)
      { atbbSouthEast(t,adj,fdist); };
    return 0;

  }

private:

  ThermalBounceback2d() : BoundaryDescriptor2d("ThermalBounceback2d") {}

};

#endif
