#include "adjointZouHe2d.hh"

/* Velocity boundaries */

void AdjointIncompressibleZouHeVelocityWest2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    adj[jj][ii][2] += 0.5*(adj[jj][ii][7] - adj[jj][ii][1]);
    adj[jj][ii][3] += adj[jj][ii][7];
    adj[jj][ii][4] += adj[jj][ii][8];
    adj[jj][ii][5] += adj[jj][ii][1];
    adj[jj][ii][6] += 0.5*(adj[jj][ii][1] - adj[jj][ii][7]);

    adj[jj][ii][1] = 0.;
    adj[jj][ii][7] = 0.;
    adj[jj][ii][8] = 0.;
  }
}

template class AdjointIncompressibleZouHeVelocityWest2d<D2Q9>;

void AdjointIncompressibleZouHeVelocityNorth2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    adj[jj][ii][4] += 0.5*(adj[jj][ii][1] - adj[jj][ii][3]);
    adj[jj][ii][5] += adj[jj][ii][1];
    adj[jj][ii][6] += adj[jj][ii][2];
    adj[jj][ii][7] += adj[jj][ii][3];
    adj[jj][ii][8] += 0.5*(adj[jj][ii][3] - adj[jj][ii][1]);

    adj[jj][ii][1] = 0.;
    adj[jj][ii][2] = 0.;
    adj[jj][ii][3] = 0.;
  }
}

template class AdjointIncompressibleZouHeVelocityNorth2d<D2Q9>;

void AdjointIncompressibleZouHeVelocityEast2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    adj[jj][ii][1] += adj[jj][ii][5];
    adj[jj][ii][2] += 0.5*(adj[jj][ii][5] - adj[jj][ii][3]);
    adj[jj][ii][6] += 0.5*(adj[jj][ii][3] - adj[jj][ii][5]);
    adj[jj][ii][7] += adj[jj][ii][3];
    adj[jj][ii][8] += adj[jj][ii][4];

    adj[jj][ii][3] = 0.;
    adj[jj][ii][4] = 0.;
    adj[jj][ii][5] = 0.;
  }
}

template class AdjointIncompressibleZouHeVelocityEast2d<D2Q9>;

/* Pressure boundaries */

void AdjointIncompressibleZouHePressureEast2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscScalar adjU;

  for (auto jj : boundingBox.yRange){
    adjU = -(1./6.)*adj[jj][ii][3] - (2./3.)*adj[jj][ii][4] -
      (1./6.)*adj[jj][ii][5];
    adj[jj][ii][0] += adjU;
    adj[jj][ii][1] += 2.*adjU + adj[jj][ii][5];
    adj[jj][ii][2] += adjU + 0.5*(adj[jj][ii][5] - adj[jj][ii][3]);
    adj[jj][ii][6] += adjU + 0.5*(adj[jj][ii][3] - adj[jj][ii][5]);
    adj[jj][ii][7] += 2.*adjU + adj[jj][ii][3];
    adj[jj][ii][8] += 2.*adjU + adj[jj][ii][4];

    adj[jj][ii][3] = 0.;
    adj[jj][ii][4] = 0.;
    adj[jj][ii][5] = 0.;
  }
}

template class AdjointIncompressibleZouHePressureEast2d<D2Q9>;

void AdjointIncompressibleZouHePressureWest2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscScalar adjU;

  for (auto jj : boundingBox.yRange){
    adjU = (1./6.)*adj[jj][ii][1] + (1./6.)*adj[jj][ii][7] + (2./3.)*adj[jj][ii][8];
    adj[jj][ii][0] -= adjU;
    adj[jj][ii][2] -= adjU - 0.5*(adj[jj][ii][7] - adj[jj][ii][1]);
    adj[jj][ii][3] -= 2.*adjU - adj[jj][ii][7];
    adj[jj][ii][4] -= 2.*adjU - adj[jj][ii][8];
    adj[jj][ii][5] -= 2.*adjU - adj[jj][ii][1];
    adj[jj][ii][6] -= adjU - 0.5*(adj[jj][ii][1] - adj[jj][ii][7]);

    adj[jj][ii][1] = 0.;
    adj[jj][ii][7] = 0.;
    adj[jj][ii][8] = 0.;
  }
}

template class AdjointIncompressibleZouHePressureWest2d<D2Q9>;

void AdjointIncompressibleZouHeNeumannEast2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt dd;
  PetscScalar *adjN1, *adjN2;

  for (auto jj : boundingBox.yRange){
    adjN1 = adj[jj][ii-1];
    adjN2 = adj[jj][ii-2];
    for (dd = 0; dd < D2Q9::numDOF; ++dd){
      adjN1[dd] -= D2Q9::ex[dd]*((2./9.)*(adj[jj][ii][3] + adj[jj][ii][5])
                                 + (8./9.)*adj[jj][ii][4]);
      adjN2[dd] += D2Q9::ex[dd]*((1./18.)*(adj[jj][ii][3] + adj[jj][ii][5])
                                 + (2./9.)*adj[jj][ii][4]);
    }
    adj[jj][ii][1] += adj[jj][ii][5];
    adj[jj][ii][2] += 0.5*(adj[jj][ii][5] - adj[jj][ii][3]);
    adj[jj][ii][6] += 0.5*(adj[jj][ii][3] - adj[jj][ii][5]);
    adj[jj][ii][7] += adj[jj][ii][3];
    adj[jj][ii][8] += adj[jj][ii][4];

    adj[jj][ii][3] = 0.;
    adj[jj][ii][4] = 0.;
    adj[jj][ii][5] = 0.;
  }
}

template class AdjointIncompressibleZouHeNeumannEast2d<D2Q9>;

void AdjointIncompressibleZouHeNeumannWest2d<D2Q9>::execute(PetscInt, void* adj_v, void*) const
{
  auto adj = static_cast<PetscScalar***>(adj_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt dd;
  PetscScalar *adjN1, *adjN2;

  for (auto jj : boundingBox.yRange){
    adjN1 = adj[jj][ii+1];
    adjN2 = adj[jj][ii+2];
    for (dd = 0; dd < D2Q9::numDOF; ++dd){
      adjN1[dd] += D2Q9::ex[dd]*((2./9.)*(adj[jj][ii][1] + adj[jj][ii][7])
                                 + (8./9.)*adj[jj][ii][8]);
      adjN2[dd] -= D2Q9::ex[dd]*((1./18.)*(adj[jj][ii][1] + adj[jj][ii][7])
                                 + (2./9.)*adj[jj][ii][8]);
    }
    adj[jj][ii][2] += 0.5*(adj[jj][ii][7] - adj[jj][ii][1]);
    adj[jj][ii][3] += adj[jj][ii][7];
    adj[jj][ii][4] += adj[jj][ii][8];
    adj[jj][ii][5] += adj[jj][ii][1];
    adj[jj][ii][6] += 0.5*(adj[jj][ii][1] - adj[jj][ii][7]);

    adj[jj][ii][1] = 0.;
    adj[jj][ii][7] = 0.;
    adj[jj][ii][8] = 0.;
  }
}

template class AdjointIncompressibleZouHeNeumannWest2d<D2Q9>;
