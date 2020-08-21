#ifndef OUTFLOW2D
#define OUTFLOW2D

#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

template <class CollisionOperator>
class OutflowEast2d : public LBMacroFunctional {

  using Lattice = typename CollisionOperator::LatticeType;
  using Macros = typename CollisionOperator::Macros ;
public:
  void evaluate(PetscInt,const PetscScalar timeConst,
                void* mac_v, void*, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar of = 0.;

    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        of -= timeConst*invOutflowLength*mac[jj][ii][1];
      }
    }
    MPI_Allreduce(&of,value,1,MPIU_SCALAR,MPI_SUM,communicator);
  }
  void adjointCollideSource(PetscInt,const PetscScalar timeConst,
                            void* adj_v, void* fdist_v, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    if (std::is_same<Macros,IncompressibleMacros2d<Lattice>>::value){
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
  OutflowEast2d(const Box2d globalOutflow, const NewLBSolver2d& solver)
  {
    communicator = solver.communicator;
    boxIntersection(solver.getLocalBoundingBox(),globalOutflow,outflow);
    PetscInt outflowLength = globalOutflow.getNx()*globalOutflow.getNy();
    invOutflowLength = (PetscScalar) (1./outflowLength);
  }

private:
  Box2d outflow;
  PetscScalar invOutflowLength;
  MPI_Comm communicator;
};

template <class CollisionOperator>
class PositiveUxFlow2d : public LBMacroFunctional {

  using Lattice = typename CollisionOperator::LatticeType;
  using Macros = typename CollisionOperator::Macros ;
public:
  void evaluate(PetscInt,const PetscScalar timeConst,
                void* mac_v, void*, PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar of = 0.;

    for (auto jj : inflow.yRange){
      for (auto ii : inflow.xRange){
        of -= timeConst*invInflowLength*mac[jj][ii][1];
      }
    }
    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        of -= timeConst*invOutflowLength*mac[jj][ii][1];
      }
    }
    MPI_Allreduce(&of,value,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
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
  PositiveUxFlow2d(const Box2d globalInflow, const Box2d globalOutflow,
                 const Box2d localBoundingBox)
  {
    boxIntersection(localBoundingBox,globalInflow,inflow);
    boxIntersection(localBoundingBox,globalOutflow,outflow);
    PetscInt inflowLength = globalInflow.getNx()*globalInflow.getNy();
    PetscInt outflowLength = globalOutflow.getNx()*globalOutflow.getNy();
    invInflowLength = (PetscScalar) (1./inflowLength);
    invOutflowLength = (PetscScalar) (1./outflowLength);
  }
private:
  Box2d inflow,outflow;
  PetscScalar invInflowLength,invOutflowLength;
};

#endif
