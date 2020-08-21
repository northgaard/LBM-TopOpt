#ifndef META
#define META

/*
  Don't ask me questions about this one, C++ metaprogramming idioms
  are weird.
*/

/*
  These macros allow introspection.

  MAKE_HAS_MEMBER_FUNCTION creates a struct that can check for the existence
  of a member function (note that there is no checking of function signature,
  merely that a function of the given name exists).

  MAKE_HAS_TYPEDEF creates a struct which can check for the existence of a nested
  typedef, i.e. a type alias made using either using = ... or typedef ... syntax.

  Example usage:
  MAKE_HAS_MEMBER_FUNCTION(doSomething,has_do_something_method);
  if (has_do_something_method<MyFancyClass>::value){
  ...
  } else {
  ...
  }
*/

#define MAKE_HAS_MEMBER_FUNCTION(func,name)                             \
  template <class T>                                                    \
  struct name                                                           \
  {                                                                     \
    template <class,class> class checker;                               \
    template <class C>                                                  \
      static std::true_type test(checker<C,decltype(&C::func)>*);       \
    template <class C>                                                  \
      static std::true_type test(checker<C,decltype(&C::template func<>)>*); \
    template <class C>                                                  \
      static std::false_type test(...);                                 \
    using type = decltype(test<T>(nullptr));                            \
    static constexpr bool value =                                       \
      std::is_same<std::true_type,decltype(test<T>(nullptr))>::value;   \
  }

#define MAKE_HAS_TYPEDEF(def,name)                                    \
  template <class T>                                                  \
  struct name                                                         \
  {                                                                   \
    template <class,class> class checker;                             \
    template <class C>                                                \
      static std::true_type test(checker<C,typename C::def>*);       \
    template <class C>                                                \
      static std::false_type test(...);                               \
    using type = decltype(test<T>(nullptr));                          \
    static constexpr bool value =                                     \
      std::is_same<std::true_type,decltype(test<T>(nullptr))>::value; \
  }

MAKE_HAS_MEMBER_FUNCTION(getSourceAdjoint,has_source_adjoint);
MAKE_HAS_MEMBER_FUNCTION(getCodiAdjoint,has_codi_adjoint);
MAKE_HAS_MEMBER_FUNCTION(diff,has_diff_function);

#endif
