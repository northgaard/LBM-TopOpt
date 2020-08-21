#ifndef THROUGHFLOW2D
#define THROUGHFLOW2D

#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

template <class CollisionOperator>
class ThroughFlow2d : public LBMacroFunctional {

  using Lattice = typename CollisionOperator::LatticeType;
  using Macros = typename CollisionOperator::Macros;
public:

  void evaluate(PetscInt, const PetscScalar timeConst,
                void* mac_v, void*, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar iif = 0.;
    PetscScalar of = 0.;
    PetscScalar inflowTot = 0., outflowTot = 0.;

    for (auto jj : inflow.yRange){
      for (auto ii : inflow.xRange){
        iif -= timeConst*invInflowLength*mac[jj][ii][1];
      }
    }
    MPI_Allreduce(&iif,&inflowTot,1,MPIU_SCALAR,MPI_SUM,communicator);

    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        of -= timeConst*invOutflowLength*mac[jj][ii][1];
      }
    }
    MPI_Allreduce(&of,&outflowTot,1,MPIU_SCALAR,MPI_SUM,communicator);

    outflowValue += outflowTot;
    *value = inflowTot + outflowTot;
  }
  void adjointCollideSource(PetscInt, const PetscScalar timeConst,
                            void* adj_v, void* fdist_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    if (std::is_same<Macros,IncompressibleMacros2d<Lattice>>::value){

      for (auto jj : inflow.yRange){
        for (auto ii : inflow.xRange){
          for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
            adj[jj][ii][dd] -= invInflowLength*timeConst*Lattice::ex[dd];
          }
        }
      }

      for (auto jj : outflow.yRange){
        for (auto ii : outflow.xRange){
          for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
            adj[jj][ii][dd] -= invOutflowLength*timeConst*Lattice::ex[dd];
          }
        }
      }

    } else if (std::is_same<Macros,StandardMacros2d<Lattice>>::value){
      auto fdist = static_cast<PetscScalar***>(fdist_v);
      PetscScalar rho,ux,uy,invRho;

      for (auto jj : inflow.yRange){
        for (auto ii : inflow.xRange){
          Macros::compute(fdist[jj][ii],rho,ux,uy);
          invRho = 1./rho;
          for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
            adj[jj][ii][dd] -= invInflowLength*timeConst*invRho*(Lattice::ex[dd] - ux);
          }
        }
      }

      for (auto jj : outflow.yRange){
        for (auto ii : outflow.xRange){
          Macros::compute(fdist[jj][ii],rho,ux,uy);
          invRho = 1./rho;
          for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
            adj[jj][ii][dd] -= invOutflowLength*timeConst*invRho*(Lattice::ex[dd] - ux);
          }
        }
      }
    }
  }
  ThroughFlow2d(const Box2d globalInflow, const Box2d globalOutflow,
                const NewLBSolver2d& solver)
  {
    communicator = solver.communicator;
    boxIntersection(solver.getLocalBoundingBox(),globalInflow,inflow);
    boxIntersection(solver.getLocalBoundingBox(),globalOutflow,outflow);
    PetscInt inflowLength = globalInflow.getNx()*globalInflow.getNy();
    invInflowLength = (PetscScalar) (1./inflowLength);
    PetscInt outflowLength = globalOutflow.getNx()*globalOutflow.getNy();
    invOutflowLength = (PetscScalar) (1./outflowLength);
  }
  PetscScalar getOutflowValue()
  {
    return outflowValue;
  }
  void reset()
  {
    outflowValue = 0.;
  }
private:
  Box2d inflow, outflow;
  PetscScalar invInflowLength, invOutflowLength;
  mutable PetscScalar outflowValue;
  MPI_Comm communicator;
};

#endif
