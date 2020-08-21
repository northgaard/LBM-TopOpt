#include "SobolGenerator.hh"
#include <string>
using std::string;
#include "SOBOL/sobol.hpp"
#include <cmath>
#include "boost/math/special_functions/erf.hpp"
#include <type_traits>

namespace {
  template <class T>
  struct dependency : public std::false_type {};
  template <class RealType>
  void SobolDispath(int,void*,RealType*)
  {
    // The dummy class forces dependency on template parameter
    static_assert(dependency<RealType>::value,"Attempting to compile SobolGenerator with incompatible floating point type, PetscReal must be float or double.\n");
  }
  template <>
  void SobolDispath(int dimNum,void* seed_v, float* array)
  {
    auto seed = static_cast<typename soboldetail::get_seed<float>::type*>(seed_v);
    i4_sobol(dimNum,seed,array);
  }
  template <>
  void SobolDispath(int dimNum,void* seed_v, double* array)
  {
    auto seed = static_cast<typename soboldetail::get_seed<double>::type*>(seed_v);
    i8_sobol(dimNum,seed,array);
  }
}

template <class RealType>
SobolGenerator_t<RealType>::SobolGenerator_t(SobolSeedType _s) : seed(_s){}

template <class RealType>
std::vector<RealType> SobolGenerator_t<RealType>::getSequence(int dimNum)
{
  std::vector<RealType> vec;
  vec.resize(dimNum);
  SobolDispath(dimNum,&seed,&vec[0]);
  return vec;
}

template <class RealType>
std::vector<RealType> SobolGenerator_t<RealType>::
getNormalDistributedSequence(int dimNum, RealType mean, RealType standardDeviation)
{
  using boost::math::erf_inv;
  std::vector<RealType> vec = getSequence(dimNum);
  for (auto& elem : vec){
    elem = mean + standardDeviation*std::sqrt(2.)*erf_inv(2.*elem - 1.);
  }
  return vec;
}

/* Explicit instantiation */
template class SobolGenerator_t<PetscReal>;
