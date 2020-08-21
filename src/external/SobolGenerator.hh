#ifndef SOBOLGENERATOR
#define SOBOLGENERATOR

#include "petscsys.h"
#include <vector>

namespace soboldetail {
  template <class RealType>
  struct get_seed {
    using type = long long int;
  };
  template <>
  struct get_seed<float> {
    using type = int;
  };
}

using SobolSeedType = typename soboldetail::get_seed<PetscReal>::type;

template <typename RealType>
class SobolGenerator_t {

public:
  SobolGenerator_t(SobolSeedType _s);
  std::vector<RealType> getSequence(int);
  std::vector<RealType> getNormalDistributedSequence(int,RealType,RealType);
private:
  SobolSeedType seed;
};

using SobolGenerator = SobolGenerator_t<PetscReal>;

#endif
