#include "symmetry2d.hh"
#include "core/macros.hh"
#include "adjointBoundaryConditions/adjointSymmetry2d.hh"

/* D2Q9 implementation */

void SymmetrySouth2d<D2Q9>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar***>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[jj][ii][5] = fdist[jj][ii][1];
    fdist[jj][ii][6] = fdist[jj][ii][2];
    fdist[jj][ii][7] = fdist[jj][ii][3];
  }
}

PetscErrorCode SymmetrySouth2d<D2Q9>::
getAdjointBoundaryLoop(AdjointLBBoundaryLoop** adjLoop) const
{
  PetscFunctionBeginUser;
  *adjLoop = new (std::nothrow) AdjointSymmetrySouth2d<D2Q9>(boundingBox);
  CHKNEWPTR(*adjLoop);
  PetscFunctionReturn(0);
}

template class SymmetrySouth2d<D2Q9>;

/* Descriptor */

PetscErrorCode Symmetry2d<D2Q9>::getSouthBoundary(const Box2d boundaryLocation,
                                                  LBBoundaryLoop** loop) const
{
  PetscFunctionBeginUser;
  *loop = new (std::nothrow) SymmetrySouth2d<D2Q9>(boundaryLocation);
  CHKNEWPTR(*loop);
  PetscFunctionReturn(0);
}
