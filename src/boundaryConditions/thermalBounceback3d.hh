#ifndef THERMALBOUNCEBACK3D
#define THERMALBOUNCEBACK3D

#include "petsc.h"
#include "boundaryConditions/LBBoundaryLoop.hh"
#include "boundaryConditions/boundaryDescriptor3d.hh"
#include "core/geometry3d.hh"
#include "core/macros.hh"
#include "latticeBoltzmann/lattices3d.hh"

/******
       Base templates
******/

/* Faces */

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthFace3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthFace3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackWestFace3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackEastFace3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackFrontFace3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackBackFace3d : public LBBoundaryLoop {};

/* Edges */

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthWestEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthEastEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthFrontEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthBackEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthWestEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthEastEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthFrontEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthBackEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackWestFrontEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackWestBackEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackEastFrontEdge3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackEastBackEdge3d : public LBBoundaryLoop {};

/* Corners */

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthFaceNorthWestCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthFaceNorthEastCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthFaceSouthWestCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackNorthFaceSouthEastCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthFaceNorthWestCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthFaceNorthEastCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthFaceSouthWestCorner3d : public LBBoundaryLoop {};

template <class IsoLattice, class ThermalLattice>
class ThermalBouncebackSouthFaceSouthEastCorner3d : public LBBoundaryLoop {};

#endif
