#ifndef COLLISIONTRAITS3D
#define COLLISIONTRAITS3D

#include "petsc.h"

template <class Operator>
struct get_num_additional_fields {
  static constexpr PetscInt value = 0;
};

template <class Operator>
struct get_lattice {
  using type = void;
};

template <>
template <class Lattice, class U>
struct get_lattice<StandardBGK3d<Lattice,U>> {
  using type = Lattice;
};

template <>
template <class Lattice, class U, class R>
struct get_lattice<StandardMRT3d<Lattice,U,R>> {
  using type = Lattice;
};

template <>
template <class Lattice, class U, class R>
struct get_lattice<IncompressibleMRT3d<Lattice,U,R>> {
  using type = Lattice;
};

#endif
