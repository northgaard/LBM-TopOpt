#include "adjointHecht3d.hh"

/******
       D3Q19 Implementation - velocity
******/

void AdjointIncompressibleHechtVelocityWest3d<D3Q19>::
execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar****>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscScalar* a;

  for (auto kk : boundingBox.zRange){
    for (auto jj : boundingBox.yRange){
      a = adj[kk][jj][ii];
      a[1] += 0.5*(a[15] - a[7]);
      a[2] += 0.5*(a[17] - a[9]);
      a[3] += a[12];
      a[4] += 0.5*(a[15] + a[17] - (a[7] + a[9]));
      a[5] += 0.5*(a[9] + a[15] - (a[7] + a[17]));
      a[6] += a[15];
      a[8] += a[17];
      a[10] += 0.5*(a[7] - a[15]);
      a[11] += 0.5*(a[9] - a[17]);
      a[13] += 0.5*(a[7] + a[9] - (a[15] + a[17]));
      a[14] += 0.5*(a[7] + a[17] - (a[9] + a[15]));
      a[16] += a[7];
      a[18] += a[9];

      a[7] = 0.;
      a[9] = 0.;
      a[12] = 0.;
      a[15] = 0.;
      a[17] = 0.;
    }
  }
}

template class AdjointIncompressibleHechtVelocityWest3d<D3Q19>;

/******
       D3Q19 Implementation - pressure
******/

void AdjointIncompressibleHechtPressureEast3d<D3Q19>::
execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar****>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscScalar adjUx;
  PetscScalar* a;

  for (auto kk : boundingBox.zRange){
    for (auto jj : boundingBox.yRange){
      a = adj[kk][jj][ii];
      adjUx = -(1./3.)*a[3] - (1./6.)*(a[6] + a[8] + a[16] + a[18]);
      a[0] += adjUx;
      a[1] += adjUx + 0.5*(a[16] - a[6]);
      a[2] += adjUx + 0.5*(a[18] - a[8]);
      a[4] += adjUx + 0.5*(a[16] + a[18] - (a[6] + a[8]));
      a[5] += adjUx + 0.5*(a[8] + a[16] - (a[6] + a[18]));
      a[7] += 2.*adjUx + a[16];
      a[9] += 2.*adjUx + a[18];
      a[10] += adjUx + 0.5*(a[6] - a[16]);
      a[11] += adjUx + 0.5*(a[8] - a[18]);
      a[12] += 2.*adjUx + a[3];
      a[13] += adjUx + 0.5*(a[6] + a[8] - (a[16] + a[18]));
      a[14] += adjUx + 0.5*(a[6] + a[18] - (a[8] + a[16]));
      a[15] += 2.*adjUx + a[6];
      a[17] += 2.*adjUx + a[8];

      a[3] = 0.;
      a[6] = 0.;
      a[8] = 0.;
      a[16] = 0.;
      a[18] = 0.;
    }
  }
}

template class AdjointIncompressibleHechtPressureEast3d<D3Q19>;
