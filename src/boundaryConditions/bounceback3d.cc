#include "bounceback3d.hh"

/******
       D3Q19 implementation
******/

/* Boundary loops -- faces */
void OnGridBouncebackNorthFace3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    for (auto ii : boundingBox.xRange){
      fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
      fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
      fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
      fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
      fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
      std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
      std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
      std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
      std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
    }
  }
}

void OnGridBouncebackSouthFace3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    for (auto ii : boundingBox.xRange){
      fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
      fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
      fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
      fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
      fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
      std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
      std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
      std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
      std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
    }
  }
}

void OnGridBouncebackWestFace3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    for (auto jj : boundingBox.yRange){
      fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
      fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
      fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
      fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
      fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
      std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
      std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
      std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
      std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
    }
  }
}

void OnGridBouncebackEastFace3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    for (auto jj : boundingBox.yRange){
      fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
      fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
      fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
      fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
      fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
      std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
      std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
      std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
      std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
    }
  }
}

void OnGridBouncebackFrontFace3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    for (auto ii : boundingBox.xRange){
      fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
      fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
      fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
      fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
      fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
      std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
      std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
      std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
      std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
    }
  }
}

void OnGridBouncebackBackFace3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    for (auto ii : boundingBox.xRange){
      fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
      fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
      fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
      fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
      fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
      std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
      std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
      std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
      std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
    }
  }
}

/* Boundary loops -- edges */
void OnGridBouncebackNorthWestEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
    fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
    fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
    fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
    fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
    fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
    fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
    std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
    std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
  }
}

void OnGridBouncebackNorthEastEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
    fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
    fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
    fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
    fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
    fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
    fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
    std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
    std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
  }
}

void OnGridBouncebackNorthFrontEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
    fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
    fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
    fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
    fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
    fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
    fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
    std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
    std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
  }
}

void OnGridBouncebackNorthBackEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
    fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
    fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
    fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
    fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
    fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
    fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
    std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
    std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
  }
}

void OnGridBouncebackSouthWestEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
    fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
    fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
    fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
    fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
    fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
    fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
    std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
    std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
  }
}

void OnGridBouncebackSouthEastEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();

  for (auto kk : boundingBox.zRange){
    fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
    fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
    fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
    fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
    fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
    fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
    fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
    std::swap(fdist[kk][jj][ii][1],fdist[kk][jj][ii][10]);
    std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
  }
}

void OnGridBouncebackSouthFrontEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
    fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
    fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
    fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
    fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
    fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
    fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
    std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
    std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
  }
}

void OnGridBouncebackSouthBackEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto ii : boundingBox.xRange){
    fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
    fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
    fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
    fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
    fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
    fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
    fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
    std::swap(fdist[kk][jj][ii][3],fdist[kk][jj][ii][12]);
    std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
  }
}

void OnGridBouncebackWestFrontEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
    fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
    fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
    fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
    fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
    fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
    fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
    std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
    std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
  }
}

void OnGridBouncebackWestBackEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
    fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
    fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
    fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
    fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
    fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
    fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
    std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
    std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
  }
}

void OnGridBouncebackEastFrontEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
    fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
    fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
    fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
    fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
    fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
    fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
    std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
    std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
  }
}

void OnGridBouncebackEastBackEdge3d<D3Q19>::execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  for (auto jj : boundingBox.yRange){
    fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
    fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
    fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
    fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
    fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
    fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
    fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
    std::swap(fdist[kk][jj][ii][2],fdist[kk][jj][ii][11]);
    std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
  }
}

/* Boundary loops -- corners */
void OnGridBouncebackNorthFaceNorthWestCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
  fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
  fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
  fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
  fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
  fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
  std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
  std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
  std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
}

void OnGridBouncebackNorthFaceNorthEastCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
  fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
  fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
  fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
  fdist[kk][jj][ii][14] = fdist[kk][jj][ii][5];
  fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
  std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
  std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
  std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
}

void OnGridBouncebackNorthFaceSouthWestCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
  fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
  fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
  fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
  fdist[kk][jj][ii][9] = fdist[kk][jj][ii][18];
  fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
  std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
  std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
  std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
}

void OnGridBouncebackNorthFaceSouthEastCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
  fdist[kk][jj][ii][2] = fdist[kk][jj][ii][11];
  fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
  fdist[kk][jj][ii][4] = fdist[kk][jj][ii][13];
  fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
  fdist[kk][jj][ii][8] = fdist[kk][jj][ii][17];
  std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
  std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
  std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
}

void OnGridBouncebackSouthFaceNorthWestCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
  fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
  fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
  fdist[kk][jj][ii][15] = fdist[kk][jj][ii][6];
  fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
  fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
  std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
  std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
  std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
}

void OnGridBouncebackSouthFaceNorthEastCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
  fdist[kk][jj][ii][10] = fdist[kk][jj][ii][1];
  fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
  fdist[kk][jj][ii][13] = fdist[kk][jj][ii][4];
  fdist[kk][jj][ii][16] = fdist[kk][jj][ii][7];
  fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
  std::swap(fdist[kk][jj][ii][5],fdist[kk][jj][ii][14]);
  std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
  std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
}

void OnGridBouncebackSouthFaceSouthWestCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
  fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
  fdist[kk][jj][ii][7] = fdist[kk][jj][ii][16];
  fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
  fdist[kk][jj][ii][12] = fdist[kk][jj][ii][3];
  fdist[kk][jj][ii][17] = fdist[kk][jj][ii][8];
  std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
  std::swap(fdist[kk][jj][ii][6],fdist[kk][jj][ii][15]);
  std::swap(fdist[kk][jj][ii][9],fdist[kk][jj][ii][18]);
}

void OnGridBouncebackSouthFaceSouthEastCorner3d<D3Q19>::
execute(PetscInt, void* fdist_v) const
{
  auto fdist = static_cast<PetscScalar****>(fdist_v);
  PetscInt ii = boundingBox.xRange.getBeginId();
  PetscInt jj = boundingBox.yRange.getBeginId();
  PetscInt kk = boundingBox.zRange.getBeginId();

  fdist[kk][jj][ii][1] = fdist[kk][jj][ii][10];
  fdist[kk][jj][ii][3] = fdist[kk][jj][ii][12];
  fdist[kk][jj][ii][5] = fdist[kk][jj][ii][14];
  fdist[kk][jj][ii][6] = fdist[kk][jj][ii][15];
  fdist[kk][jj][ii][11] = fdist[kk][jj][ii][2];
  fdist[kk][jj][ii][18] = fdist[kk][jj][ii][9];
  std::swap(fdist[kk][jj][ii][4],fdist[kk][jj][ii][13]);
  std::swap(fdist[kk][jj][ii][7],fdist[kk][jj][ii][16]);
  std::swap(fdist[kk][jj][ii][8],fdist[kk][jj][ii][17]);
}

/******
       Explicit instantiation
******/

template class OnGridBouncebackNorthFace3d<D3Q19>;
template class OnGridBouncebackSouthFace3d<D3Q19>;
template class OnGridBouncebackWestFace3d<D3Q19>;
template class OnGridBouncebackEastFace3d<D3Q19>;
template class OnGridBouncebackFrontFace3d<D3Q19>;
template class OnGridBouncebackBackFace3d<D3Q19>;

template class OnGridBouncebackNorthWestEdge3d<D3Q19>;
template class OnGridBouncebackNorthEastEdge3d<D3Q19>;
template class OnGridBouncebackNorthFrontEdge3d<D3Q19>;
template class OnGridBouncebackNorthBackEdge3d<D3Q19>;
template class OnGridBouncebackSouthWestEdge3d<D3Q19>;
template class OnGridBouncebackSouthEastEdge3d<D3Q19>;
template class OnGridBouncebackSouthFrontEdge3d<D3Q19>;
template class OnGridBouncebackSouthBackEdge3d<D3Q19>;
template class OnGridBouncebackWestFrontEdge3d<D3Q19>;
template class OnGridBouncebackWestBackEdge3d<D3Q19>;
template class OnGridBouncebackEastFrontEdge3d<D3Q19>;
template class OnGridBouncebackEastBackEdge3d<D3Q19>;

template class OnGridBouncebackNorthFaceNorthWestCorner3d<D3Q19>;
template class OnGridBouncebackNorthFaceNorthEastCorner3d<D3Q19>;
template class OnGridBouncebackNorthFaceSouthWestCorner3d<D3Q19>;
template class OnGridBouncebackNorthFaceSouthEastCorner3d<D3Q19>;
template class OnGridBouncebackSouthFaceNorthWestCorner3d<D3Q19>;
template class OnGridBouncebackSouthFaceNorthEastCorner3d<D3Q19>;
template class OnGridBouncebackSouthFaceSouthWestCorner3d<D3Q19>;
template class OnGridBouncebackSouthFaceSouthEastCorner3d<D3Q19>;
