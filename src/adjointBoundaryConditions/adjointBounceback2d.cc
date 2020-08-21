#include "adjointBounceback2d.hh"
#include <utility>

/* D2Q9 implementation */

void AdjointOnGridBouncebackNorth2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    adj[jj][ii][5] += adj[jj][ii][1];
    adj[jj][ii][6] += adj[jj][ii][2];
    adj[jj][ii][7] += adj[jj][ii][3];
    std::swap(adj[jj][ii][4],adj[jj][ii][8]);

    adj[jj][ii][1] = 0.;
    adj[jj][ii][2] = 0.;
    adj[jj][ii][3] = 0.;
  }
}

void AdjointOnGridBouncebackSouth2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    adj[jj][ii][1] += adj[jj][ii][5];
    adj[jj][ii][2] += adj[jj][ii][6];
    adj[jj][ii][3] += adj[jj][ii][7];
    std::swap(adj[jj][ii][4],adj[jj][ii][8]);

    adj[jj][ii][5] = 0.;
    adj[jj][ii][6] = 0.;
    adj[jj][ii][7] = 0.;
  }
}

void AdjointOnGridBouncebackWest2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    adj[jj][ii][3] += adj[jj][ii][7];
    adj[jj][ii][4] += adj[jj][ii][8];
    adj[jj][ii][5] += adj[jj][ii][1];
    std::swap(adj[jj][ii][2],adj[jj][ii][6]);

    adj[jj][ii][1] = 0.;
    adj[jj][ii][7] = 0.;
    adj[jj][ii][8] = 0.;
  }
}

void AdjointOnGridBouncebackEast2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    adj[jj][ii][1] += adj[jj][ii][5];
    adj[jj][ii][7] += adj[jj][ii][3];
    adj[jj][ii][8] += adj[jj][ii][4];
    std::swap(adj[jj][ii][2],adj[jj][ii][6]);

    adj[jj][ii][3] = 0.;
    adj[jj][ii][4] = 0.;
    adj[jj][ii][5] = 0.;
  }
}

void AdjointOnGridBouncebackNorthWest2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  adj[jj][ii][4] += adj[jj][ii][8];
  adj[jj][ii][5] += adj[jj][ii][1];
  adj[jj][ii][6] += adj[jj][ii][2];
  std::swap(adj[jj][ii][3],adj[jj][ii][7]);

  adj[jj][ii][1] = 0.;
  adj[jj][ii][2] = 0.;
  adj[jj][ii][8] = 0.;
}

void AdjointOnGridBouncebackNorthEast2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  adj[jj][ii][6] += adj[jj][ii][2];
  adj[jj][ii][7] += adj[jj][ii][3];
  adj[jj][ii][8] += adj[jj][ii][4];
  std::swap(adj[jj][ii][1],adj[jj][ii][5]);

  adj[jj][ii][2] = 0.;
  adj[jj][ii][3] = 0.;
  adj[jj][ii][4] = 0.;
}

void AdjointOnGridBouncebackSouthWest2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  adj[jj][ii][2] += adj[jj][ii][6];
  adj[jj][ii][3] += adj[jj][ii][7];
  adj[jj][ii][4] += adj[jj][ii][8];
  std::swap(adj[jj][ii][1],adj[jj][ii][5]);

  adj[jj][ii][6] = 0.;
  adj[jj][ii][7] = 0.;
  adj[jj][ii][8] = 0.;
}

void AdjointOnGridBouncebackSouthEast2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  PetscScalar*** adj = (PetscScalar***) adj_v;
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  adj[jj][ii][1] += adj[jj][ii][5];
  adj[jj][ii][2] += adj[jj][ii][6];
  adj[jj][ii][8] += adj[jj][ii][4];
  std::swap(adj[jj][ii][3],adj[jj][ii][7]);

  adj[jj][ii][4] = 0.;
  adj[jj][ii][5] = 0.;
  adj[jj][ii][6] = 0.;
}
