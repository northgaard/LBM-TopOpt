#ifndef CORETRAITS
#define CORETRAITS

#include "petsc.h"
#include <type_traits>

template <PetscInt N>
struct integer_to_bool_type {
  using type = std::true_type;
};

template <>
struct integer_to_bool_type<0> {
  using type = std::false_type;
};

#endif
