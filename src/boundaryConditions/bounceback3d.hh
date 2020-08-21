#ifndef BOUNCEBACK3D
#define BOUNCEBACK3D

#include "petsc.h"
#include "LBBoundaryLoop.hh"
#include "core/geometry3d.hh"
#include "core/macros.hh"
#include "latticeBoltzmann/lattices3d.hh"
#include "boundaryConditions/boundaryDescriptor3d.hh"
#include <utility>

/******
       Base templates
******/

/* Faces */

template <class Lattice>
class OnGridBouncebackNorthFace3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthFace3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackWestFace3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackEastFace3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackFrontFace3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackBackFace3d : public LBBoundaryLoop {};

/* Edges */

template <class Lattice>
class OnGridBouncebackNorthWestEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackNorthEastEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackNorthFrontEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackNorthBackEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthWestEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthEastEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthFrontEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthBackEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackWestFrontEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackWestBackEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackEastFrontEdge3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackEastBackEdge3d : public LBBoundaryLoop {};

/* Corners */

template <class Lattice>
class OnGridBouncebackNorthFaceNorthWestCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackNorthFaceNorthEastCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackNorthFaceSouthWestCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackNorthFaceSouthEastCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthFaceNorthWestCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthFaceNorthEastCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthFaceSouthWestCorner3d : public LBBoundaryLoop {};

template <class Lattice>
class OnGridBouncebackSouthFaceSouthEastCorner3d : public LBBoundaryLoop {};

/******
       D3Q19 implementation
******/

template <>
class OnGridBouncebackNorthFace3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthFace3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                           "NorthFace"), boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthFace3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthFace3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                           "SouthFace"), boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackWestFace3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackWestFace3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                           "WestFace"), boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackEastFace3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackEastFace3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                           "EastFace"), boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackFrontFace3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackFrontFace3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                           "FrontFace"), boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackBackFace3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackBackFace3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                           "BackFace"), boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthWestEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthWestEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "NorthWestEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthEastEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthEastEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "NorthEastEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthFrontEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthFrontEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                                "NorthFrontEdge"),
                                                 boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthBackEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthBackEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "NorthBackEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthWestEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthWestEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "SouthWestEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthEastEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthEastEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "SouthEastEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthFrontEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthFrontEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                                "SouthFrontEdge"),
                                                 boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthBackEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthBackEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                                "SouthBackEdge"),
                                                 boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackWestFrontEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackWestFrontEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "WestFrontEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackWestBackEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackWestBackEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                              "WestBackEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackEastFrontEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackEastFrontEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                              "EastFrontEdge"),
                                               boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackEastBackEdge3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackEastBackEdge3d(Box3d _box) : LBBoundaryLoop("OnGridBounceback3d",
                                                               "EastBackEdge"),
                                                boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthFaceNorthWestCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthFaceNorthWestCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","NorthFaceNorthWestCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthFaceNorthEastCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthFaceNorthEastCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","NorthFaceNorthEastCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthFaceSouthWestCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthFaceSouthWestCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","NorthFaceSouthWestCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackNorthFaceSouthEastCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackNorthFaceSouthEastCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","NorthFaceSouthEastCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthFaceNorthWestCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthFaceNorthWestCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","SouthFaceNorthWestCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthFaceNorthEastCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthFaceNorthEastCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","SouthFaceNorthEastCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthFaceSouthWestCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthFaceSouthWestCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","SouthFaceSouthWestCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

template <>
class OnGridBouncebackSouthFaceSouthEastCorner3d<D3Q19> : public LBBoundaryLoop {

public:
  OnGridBouncebackSouthFaceSouthEastCorner3d(Box3d _box) :
    LBBoundaryLoop("OnGridBounceback3d","SouthFaceSouthEastCorner"),
    boundingBox(_box){}
  void execute(PetscInt,void*) const override;
private:
  Box3d boundingBox;
};

/******
       Boundary descriptor
******/

template <class Lattice>
class OnGridBounceback3d : public BoundaryDescriptor3d {

public:
  OnGridBounceback3d() : BoundaryDescriptor3d("OnGridBounceback3d"){}
  PetscErrorCode getNorthBoundary(const Box3d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthFace3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthBoundary(const Box3d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthFace3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestBoundary(const Box3d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackWestFace3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBoundary(const Box3d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackEastFace3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getFrontBoundary(const Box3d boundaryLocation,
                                  LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackFrontFace3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getBackBoundary(const Box3d boundaryLocation,
                                 LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackBackFace3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthWestEdgeBoundary(const Box3d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthWestEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthEastEdgeBoundary(const Box3d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthEastEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthFrontEdgeBoundary(const Box3d boundaryLocation,
                                       LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthFrontEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthBackEdgeBoundary(const Box3d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthBackEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthWestEdgeBoundary(const Box3d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthWestEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthEastEdgeBoundary(const Box3d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthEastEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthFrontEdgeBoundary(const Box3d boundaryLocation,
                                       LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthFrontEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthBackEdgeBoundary(const Box3d boundaryLocation,
                                      LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthBackEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestFrontEdgeBoundary(const Box3d boundaryLocation,
                                          LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackWestFrontEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getWestBackEdgeBoundary(const Box3d boundaryLocation,
                                         LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackWestBackEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastFrontEdgeBoundary(const Box3d boundaryLocation,
                                          LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackEastFrontEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getEastBackEdgeBoundary(const Box3d boundaryLocation,
                                         LBBoundaryLoop** loop) const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackEastBackEdge3d<Lattice>(boundaryLocation);
    CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthFaceNorthWestCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthFaceNorthWestCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthFaceNorthEastCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthFaceNorthEastCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthFaceSouthWestCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthFaceSouthWestCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getNorthFaceSouthEastCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackNorthFaceSouthEastCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthFaceNorthWestCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthFaceNorthWestCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthFaceNorthEastCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthFaceNorthEastCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthFaceSouthWestCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthFaceSouthWestCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
  PetscErrorCode getSouthFaceSouthEastCornerBoundary(const Box3d boundaryLocation,
                                                     LBBoundaryLoop** loop)
    const override
  {
    PetscFunctionBeginUser;
    *loop = new (std::nothrow) OnGridBouncebackSouthFaceSouthEastCorner3d<Lattice>
      (boundaryLocation); CHKNEWPTR(*loop);
    PetscFunctionReturn(0);
  }
};

#endif
