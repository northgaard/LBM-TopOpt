#ifndef PRESSUREDROP2D
#define PRESSUREDROP2D

#include "petsc.h"
#include "functionals/LBMacroFunctional.hh"
#include "core/geometry2d.hh"
#include "latticeBoltzmann/NewLBSolver2d.hh"

template <class CollisionOperator>
class PressureDrop2d : public LBMacroFunctional {

  using Lattice = typename CollisionOperator::LatticeType;
public:
  void evaluate(PetscInt,const PetscScalar timeConst,
                void* mac_v,void*,PetscScalar* value) const override
  {
    auto mac = static_cast<PetscScalar***>(mac_v);
    *value = 0.;
    PetscScalar ps = 0.;

    for (auto jj : inflow.yRange){
      for (auto ii : inflow.xRange){
        ps += timeConst*invInflowLength*mac[jj][ii][0];
      }
    }
    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        ps -= timeConst*invOutflowLength*mac[jj][ii][0];
      }
    }
    MPI_Allreduce(&ps,value,1,MPIU_SCALAR,MPI_SUM,communicator);
  }
  void adjointCollideSource(PetscInt,const PetscScalar timeConst,
                            void* adj_v, void*, void*) const override
  {
    auto adj = static_cast<PetscScalar***>(adj_v);
    for (auto jj : inflow.yRange){
      for (auto ii : inflow.xRange){
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          adj[jj][ii][dd] += invInflowLength*timeConst;
        }
      }
    }

    for (auto jj : outflow.yRange){
      for (auto ii : outflow.xRange){
        for (size_t dd = 0; dd < Lattice::numDOF; ++dd){
          adj[jj][ii][dd] -= invOutflowLength*timeConst;
        }
      }
    }
  }
  PressureDrop2d(const Box2d globalInflow, const Box2d globalOutflow,
                 const NewLBSolver2d& solver)
  {
    communicator = solver.communicator;
    boxIntersection(solver.getLocalBoundingBox(),globalInflow,inflow);
    boxIntersection(solver.getLocalBoundingBox(),globalOutflow,outflow);
    PetscInt inflowLength = globalInflow.getNx()*globalInflow.getNy();
    PetscInt outflowLength = globalOutflow.getNx()*globalOutflow.getNy();
    invInflowLength = (PetscScalar) (1./inflowLength);
    invOutflowLength = (PetscScalar) (1./outflowLength);
  }

private:
  Box2d inflow;
  Box2d outflow;
  PetscScalar invInflowLength,invOutflowLength;
  MPI_Comm communicator;
};

#endif
