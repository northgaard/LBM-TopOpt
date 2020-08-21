#ifndef COLLISIONTRAITS2D
#define COLLISIONTRAITS2D

#include "petsc.h"

template <class Operator>
struct has_static_equilibrium {
  static constexpr bool value = true;
};

template <>
template <class L, class U>
struct has_static_equilibrium<IncompressibleBGKUnit2d<L,U>> {
  static constexpr bool value = false;
};

template <>
template <class L, class U>
struct has_static_equilibrium<StandardBGKUnit2d<L,U>> {
  static constexpr bool value = false;
};

template <>
template <class BaseOperator, class I>
struct has_static_equilibrium<PartialBouncebackCollision<BaseOperator,I>> {
  static constexpr bool value = has_static_equilibrium<BaseOperator>::value;
};

#endif
