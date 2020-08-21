#ifndef LATTICES2D
#define LATTICES2D

/* Specific lattices */

/* IMPORTANT: The ordering of the lattice velocities is NOT irrelevant.
   It needs to be in a certain order to enable correctly performing
   collision and streaming in one loop.

   See: Jonas Latt -- "How to implement you DdQq dynamics with only
   q variables per node (instead of 2q)"
*/

#include "petsc.h"

struct D2Q9 {
  constexpr static PetscInt numDOF = 9;
  constexpr static PetscInt half = (numDOF - 1)/2;
  constexpr static PetscInt numMacros = 3;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar cs = 0.5773502691896258; // 1.0/sqrt(3.0)
  constexpr static PetscScalar csSq = 1./3.;
  constexpr static PetscScalar csInv = 1.7320508075688772; // sqrt(3.0)
  constexpr static PetscScalar csSqInv = 3.;

  constexpr static PetscInt ex[numDOF] = {0,1,0,-1,-1,-1,0,1,1};
  constexpr static PetscInt ey[numDOF] = {0,-1,-1,-1,0,1,1,1,0};
  constexpr static PetscInt opp[numDOF] = {0,5,6,7,8,1,2,3,4};
  constexpr static PetscScalar weights[numDOF] = {4./9.,1./36.,1./9.,1./36.,1./9.,
					    1./36.,1./9.,1./36.,1./9.};
};

struct D2Q4 {
  constexpr static PetscInt numDOF = 4;
  constexpr static PetscInt half = numDOF / 2;
  constexpr static PetscInt numMacros = 1;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar cs = 0.7071067811865475; // 1.0/sqrt(2.0)
  constexpr static PetscScalar csSq = 1./2.;
  constexpr static PetscScalar csInv = 1.4142135623730951; // sqrt(2.0)
  constexpr static PetscScalar csSqInv = 2.;

  constexpr static PetscInt ex[numDOF] = {0,-1,0,1};
  constexpr static PetscInt ey[numDOF] = {-1,0,1,0};
  constexpr static PetscInt opp[numDOF] = {2,3,0,1};
  constexpr static PetscScalar weights[numDOF] = {1./4.,1./4.,1./4.,1./4.};
};

struct D2Q5 {
  constexpr static PetscInt numDOF = 5;
  constexpr static PetscInt half = (numDOF - 1) / 2;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar cs = 0.5773502691896258; // 1.0/sqrt(3.0)
  constexpr static PetscScalar csSq = 1./3.;
  constexpr static PetscScalar csInv = 1.7320508075688772; // sqrt(3.0)
  constexpr static PetscScalar csSqInv = 3.;

  constexpr static PetscInt ex[numDOF] = {0,0,-1,0,1};
  constexpr static PetscInt ey[numDOF] = {0,-1,0,1,0};
  constexpr static PetscInt opp[numDOF] = {0,3,4,1,2};
  constexpr static PetscScalar weights[numDOF] = {1./3.,1./6.,1./6.,1./6.,1./6.};
};

/*
  For some reason the double thermal MRT uses a different value of the speed of sound
  from the BGK (at least that's how it seems from the papers). I don't understand the theory
  well enough to explain why this is. For now it gets its own lattice like the special little
  snowflake it is.
*/

struct D2Q5_TMRT {
  constexpr static PetscInt numDOF = 5;
  constexpr static PetscInt half = (numDOF - 1) / 2;
  constexpr static PetscInt connectivity = 1;
  constexpr static PetscScalar cs = 0.4472135954999579; // 1.0/sqrt(5.0)
  constexpr static PetscScalar csSq = 1./5.;
  constexpr static PetscScalar csInv = 2.23606797749979; // sqrt(5.0)
  constexpr static PetscScalar csSqInv = 5.;

  constexpr static PetscInt ex[numDOF] = {0,0,-1,0,1};
  constexpr static PetscInt ey[numDOF] = {0,-1,0,1,0};
  constexpr static PetscInt opp[numDOF] = {0,3,4,1,2};
  constexpr static PetscScalar weights[numDOF] = {1./3.,1./6.,1./6.,1./6.,1./6.};
};

#endif
