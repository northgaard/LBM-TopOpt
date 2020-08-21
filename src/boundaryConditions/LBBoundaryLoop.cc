#include "LBBoundaryLoop.hh"

LBBoundaryLoop::LBBoundaryLoop(const std::string& _bn, const std::string& _ori)
  : boundaryName(_bn), orientation(_ori)
{}

LBBoundaryLoop::LBBoundaryLoop(std::string&& _bn, std::string&& _ori)
  : boundaryName(std::move(_bn)), orientation(std::move(_ori))
{}

PetscErrorCode LBBoundaryLoop::getAdjointBoundaryLoop(AdjointLBBoundaryLoop**) const
{
  PetscFunctionBeginUser;
  SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,
           "Adjoint boundary is not implemented for the %s orientation of the %s boundary type\n",
           orientation.c_str(),boundaryName.c_str());
  PetscFunctionReturn(0);
}
