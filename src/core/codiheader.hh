#ifndef CODIHEADER
#define CODIHEADER

#include "codi.hpp"
#ifdef TOPLB_HAS_ADEPT
#include "adept.hpp"
#ifndef ADEPT_STACK_THREAD_UNSAFE
#define ADEPT_STACK_THREAD_UNSAFE
#endif
#endif

template <class Real>
using CodiReverseType = codi::ActiveReal<codi::JacobiTape<codi::ChunkTapeTypes<
                                                            Real,
                                                            codi::LinearIndexHandler<PetscInt>>>>;
using ReversePetscScalar = CodiReverseType<PetscScalar>;

#endif
