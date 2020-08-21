#include "bounceback2d.hh"
#include "core/macros.hh"
#include "adjointBoundaryConditions/adjointBounceback2d.hh"

/* D2Q9 implementation */

void OnGridBouncebackNorth2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[jj][ii][1] = fdist[jj][ii][5];
    fdist[jj][ii][2] = fdist[jj][ii][6];
    fdist[jj][ii][3] = fdist[jj][ii][7];
    std::swap(fdist[jj][ii][4],fdist[jj][ii][8]);
  }
}

PetscErrorCode OnGridBouncebackNorth2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackNorth2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackNorth2d<D2Q9>;

void OnGridBouncebackSouth2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[jj][ii][5] = fdist[jj][ii][1];
    fdist[jj][ii][6] = fdist[jj][ii][2];
    fdist[jj][ii][7] = fdist[jj][ii][3];
    std::swap(fdist[jj][ii][4],fdist[jj][ii][8]);
  }
}

PetscErrorCode OnGridBouncebackSouth2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackSouth2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackSouth2d<D2Q9>;

void OnGridBouncebackWest2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    fdist[jj][ii][1] = fdist[jj][ii][5];
    fdist[jj][ii][7] = fdist[jj][ii][3];
    fdist[jj][ii][8] = fdist[jj][ii][4];
    std::swap(fdist[jj][ii][2],fdist[jj][ii][6]);
  }
}

PetscErrorCode OnGridBouncebackWest2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackWest2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackWest2d<D2Q9>;

void OnGridBouncebackEast2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    fdist[jj][ii][3] = fdist[jj][ii][7];
    fdist[jj][ii][4] = fdist[jj][ii][8];
    fdist[jj][ii][5] = fdist[jj][ii][1];
    std::swap(fdist[jj][ii][2],fdist[jj][ii][6]);
  }
}

PetscErrorCode OnGridBouncebackEast2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackEast2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackEast2d<D2Q9>;

void OnGridBouncebackNorthWest2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  fdist[jj][ii][1] = fdist[jj][ii][5];
  fdist[jj][ii][2] = fdist[jj][ii][6];
  fdist[jj][ii][8] = fdist[jj][ii][4];
  std::swap(fdist[jj][ii][3],fdist[jj][ii][7]);
}

PetscErrorCode OnGridBouncebackNorthWest2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackNorthWest2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackNorthWest2d<D2Q9>;

void OnGridBouncebackNorthEast2d<D2Q9>::execute(PetscInt timestep, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  fdist[jj][ii][2] = fdist[jj][ii][6];
  fdist[jj][ii][3] = fdist[jj][ii][7];
  fdist[jj][ii][4] = fdist[jj][ii][8];
  std::swap(fdist[jj][ii][1],fdist[jj][ii][5]);
}

PetscErrorCode OnGridBouncebackNorthEast2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackNorthEast2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackNorthEast2d<D2Q9>;

void OnGridBouncebackSouthWest2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  fdist[jj][ii][6] = fdist[jj][ii][2];
  fdist[jj][ii][7] = fdist[jj][ii][3];
  fdist[jj][ii][8] = fdist[jj][ii][4];
  std::swap(fdist[jj][ii][1],fdist[jj][ii][5]);
}

PetscErrorCode OnGridBouncebackSouthWest2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackSouthWest2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackSouthWest2d<D2Q9>;

void OnGridBouncebackSouthEast2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  PetscScalar*** fdist = (PetscScalar***) fdist_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  fdist[jj][ii][4] = fdist[jj][ii][8];
  fdist[jj][ii][5] = fdist[jj][ii][1];
  fdist[jj][ii][6] = fdist[jj][ii][2];
  std::swap(fdist[jj][ii][3],fdist[jj][ii][7]);
}

PetscErrorCode OnGridBouncebackSouthEast2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointOnGridBouncebackSouthEast2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class OnGridBouncebackSouthEast2d<D2Q9>;

/* Descriptor */

PetscErrorCode OnGridBounceback2d<D2Q9>::getNorthBoundary(const Box2d boundaryLocation,
                                                LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackNorth2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getSouthBoundary(const Box2d boundaryLocation,
                                                          LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackSouth2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getWestBoundary(const Box2d boundaryLocation,
                                                         LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackWest2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getEastBoundary(const Box2d boundaryLocation,
                                                         LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackEast2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getNorthWestBoundary(const Box2d boundaryLocation,
                                                              LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackNorthWest2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getNorthEastBoundary(const Box2d boundaryLocation,
                                                              LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackNorthEast2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getSouthWestBoundary(const Box2d boundaryLocation,
                                                              LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackSouthWest2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

PetscErrorCode OnGridBounceback2d<D2Q9>::getSouthEastBoundary(const Box2d boundaryLocation,
                                                              LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) OnGridBouncebackSouthEast2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}

template class OnGridBounceback2d<D2Q9>;
