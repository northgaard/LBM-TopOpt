#ifndef LBBOUNDARYLOOP
#define LBBOUNDARYLOOP

#include "petsc.h"
#include "adjointBoundaryConditions/adjointLBBoundaryLoop.hh"
#include <string>

class LBBoundaryLoop {

public:
  virtual void execute(PetscInt,void*,void*) const = 0;
  virtual PetscErrorCode getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const;
  LBBoundaryLoop(const std::string&, const std::string&);
  LBBoundaryLoop(std::string&&, std::string&&);
  virtual ~LBBoundaryLoop(){}
private:
  const std::string boundaryName;
  const std::string orientation;
};

#endif
