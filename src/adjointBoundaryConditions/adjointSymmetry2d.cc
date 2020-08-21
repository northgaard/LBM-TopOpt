#include "adjointSymmetry2d.hh"

/* D2Q9 implementation */

void AdjointSymmetrySouth2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    adj[jj][ii][1] += adj[jj][ii][5];
    adj[jj][ii][2] += adj[jj][ii][6];
    adj[jj][ii][3] += adj[jj][ii][7];

    adj[jj][ii][5] = 0.;
    adj[jj][ii][6] = 0.;
    adj[jj][ii][7] = 0.;
  }
}

template class AdjointSymmetrySouth2d<D2Q9>;
