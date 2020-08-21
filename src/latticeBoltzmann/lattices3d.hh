#ifndef LATTICES3D
#define LATTICES3D

/* Specific lattices */

/* IMPORTANT: The ordering of the lattice velocities is NOT irrelevant.
   It needs to be in a certain order to enable correctly performing
   collision and streaming in one loop.

   See: Jonas Latt -- "How to implement you DdQq dynamics with only
   q variables per node (instead of 2q)"
*/

#include "petsc.h"

// TODO: Test lattice implementation.
struct D3Q15 {

  constexpr static PetscInt numDOF = 15;
  constexpr static PetscInt half = (numDOF - 1)/2;
  constexpr static PetscInt numMacros = 3;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar csSq = 1./3.;

  constexpr static PetscInt ex[numDOF] = {0,0,0,-1,-1,1,-1,1,0,0,1,1,-1,1,-1};
  constexpr static PetscInt ey[numDOF] = {0,0,-1,0,-1,-1,1,1,0,1,0,1,1,-1,-1};
  constexpr static PetscInt ez[numDOF] = {0,-1,0,0,-1,-1,-1,-1,1,0,0,1,1,1,1};
  constexpr static PetscInt opp[numDOF] = {0,8,9,10,11,12,13,14,1,2,3,4,5,6,7};
  constexpr static PetscScalar weights[numDOF] = {2./9.,1./9.,1./9.,1./9.,1./72.,1./72.,
                                                  1./72.,1./72.,1./9.,1./9.,1./9.,
                                                  1./72.,1./72.,1./72.,1./72.};
};

struct D3Q19 {

  constexpr static PetscInt numDOF = 19;
  constexpr static PetscInt half = (numDOF - 1)/2;
  constexpr static PetscInt numMacros = 3;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar csSq = 1./3.;

  constexpr static PetscInt ex[numDOF] = {0,0,0,-1,0,0,-1,1,-1,1,0,0,1,0,0,1,-1,1,-1};
  constexpr static PetscInt ey[numDOF] = {0,0,-1,0,-1,1,0,0,-1,-1,0,1,0,1,-1,0,0,1,1};
  constexpr static PetscInt ez[numDOF] = {0,-1,0,0,-1,-1,-1,-1,0,0,1,0,0,1,1,1,1,0,0};
  constexpr static PetscInt opp[numDOF] = {0,10,11,12,13,14,15,16,17,18,
                                           1,2,3,4,5,6,7,8,9};
  constexpr static PetscScalar weights[numDOF] = {1./3.,1./18.,1./18.,1./18.,
                                                  1./36.,1./36.,1./36.,
                                                  1./36.,1./36.,1./36.,
                                                  1./18.,1./18.,1./18.,
                                                  1./36.,1./36.,1./36.,
                                                  1./36.,1./36.,1./36.};
};

// TODO: Test lattice implementation
struct D3Q27 {

  constexpr static PetscInt numDOF = 27;
  constexpr static PetscInt half = (numDOF - 1)/2;
  constexpr static PetscInt numMacros = 3;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar csSq = 1./3.;

  constexpr static PetscInt ex[numDOF] = {0,0,0,-1,0,0,-1,1,-1,1,-1,1,-1,1,0,0,1,
                                          0,0,1,-1,1,-1,1,-1,1,-1};
  constexpr static PetscInt ey[numDOF] = {0,0,-1,0,-1,1,0,0,-1,-1,-1,-1,1,1,0,1,0,
                                          1,-1,0,0,1,1,1,1,-1,-1};
  constexpr static PetscInt ez[numDOF] = {0,-1,0,0,-1,-1,-1,-1,0,0,-1,-1,-1,-1,1,
                                          0,0,1,1,1,1,0,0,1,1,1,1};
  constexpr static PetscInt opp[numDOF] = {0,14,15,16,17,18,19,20,21,22,23,24,25,26,
                                           1,2,3,4,5,6,7,8,9,10,11,12,13};
  constexpr static PetscScalar weights[numDOF] = {8./27.,2./27.,2./27.,2./27.,
                                                  1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
                                                  1./216.,1./216.,1./216.,1./216.,
                                                  2./27.,2./27.,2./27.,
                                                  1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
                                                  1./216.,1./216.,1./216.,1./216.};
};

#endif
