#include "finiteDifferenceCheck2d.hh"

PetscErrorCode FiniteDifferenceCheck2d::checkDomain(
    NewLBSolver2d &solver, AdjointLBSolver &adjSolver,
    const LBMacroFunctional &obj, Box2d domain, Filter *filter,
    Vec inputDesignVec) {
  PetscErrorCode ierr;
  Vec designVec, filteredVec, sensitivityVec;
  DM designGrid;
  PetscScalar baseObj, perturbedObj, fdResult, obstOld = 0.;
  PetscInt rank;
  PetscBool intersect;
  Box2d localBoundingBox = solver.getLocalBoundingBox();
  Box2d point;
  PetscScalar **obst, **sens;
  PetscFunctionBeginUser;
  ierr = solver.getDesignGrid(&designGrid);
  CHKERRQ(ierr);
  if (inputDesignVec) {
    designVec = inputDesignVec;
  } else {
    ierr = solver.createDesignVec(&designVec);
    CHKERRQ(ierr);
    ierr = VecSetRandom(designVec, nullptr);
    CHKERRQ(ierr);
    ierr = VecAbs(designVec);
    CHKERRQ(ierr);
    PetscScalar max;
    ierr = VecMax(designVec, nullptr, &max);
    CHKERRQ(ierr);
    ierr = VecScale(designVec, 1. / max);
    CHKERRQ(ierr);
  }
  if (filter) {
    ierr = VecDuplicate(designVec, &filteredVec);
    CHKERRQ(ierr);
    ierr = filter->filterDesign(designVec, filteredVec);
    CHKERRQ(ierr);
  } else {
    filteredVec = designVec;
  }
  ierr = solver.setFieldFromDesignVec(filteredVec);
  CHKERRQ(ierr);
  double t1, t2;
  t1 = MPI_Wtime();
  ierr =
      computeSensitivities(solver, adjSolver, obj, &baseObj, &sensitivityVec);
  CHKERRQ(ierr);
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "Adjoint time: %f\n", t2 - t1);
  if (filter) {
    ierr = filter->filterSensitivities(sensitivityVec, nullptr, 0);
    CHKERRQ(ierr);
  }

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  for (auto jj : domain.yRange) {
    for (auto ii : domain.xRange) {
      point = Box2d(ii, ii, jj, jj);
      intersect = doBoxesIntersect(localBoundingBox, point);

      if (intersect) {
        // Perturb design
        ierr = DMDAVecGetArray(designGrid, designVec, &obst);
        CHKERRQ(ierr);
        obstOld = obst[jj][ii];
        if (obst[jj][ii] <= 0.5) {
          obst[jj][ii] += delta;
        } else {
          obst[jj][ii] -= delta;
        }
        ierr = DMDAVecRestoreArray(designGrid, designVec, &obst);
        CHKERRQ(ierr);
      }
      // Compute perturbed objective
      if (filter) {
        filter->filterDesign(designVec, filteredVec);
        CHKERRQ(ierr);
      }
      ierr = solver.setFieldFromDesignVec(filteredVec);
      CHKERRQ(ierr);
      ierr = computeObjective(solver, obj, &perturbedObj);
      CHKERRQ(ierr);

      if (intersect) {
        // Compute finite difference sensitivity in the point
        if (obstOld <= 0.5) {
          fdResult = (perturbedObj - baseObj) / delta;
        } else {
          fdResult = (baseObj - perturbedObj) / delta;
        }
        // Restore to undisturbed state
        ierr = DMDAVecGetArray(designGrid, designVec, &obst);
        CHKERRQ(ierr);
        obst[jj][ii] = obstOld;
        ierr = DMDAVecRestoreArray(designGrid, designVec, &obst);
        CHKERRQ(ierr);
        // Compare with adjoint result
        ierr =
            DMDAVecGetArray(adjSolver.sensitivityGrid, sensitivityVec, &sens);
        CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_SELF, "--- From process %d ---\n", rank);
        PetscPrintf(PETSC_COMM_SELF, "At point: (%d,%d)\n", ii, jj);
        PetscPrintf(PETSC_COMM_SELF, "FD: %e\tAdjoint: %e\n", fdResult,
                    sens[jj][ii]);
        PetscPrintf(PETSC_COMM_SELF, "Relative: %e\n\n",
                    fdResult / sens[jj][ii]);
        ierr = DMDAVecRestoreArray(adjSolver.sensitivityGrid, sensitivityVec,
                                   &sens);
        CHKERRQ(ierr);
      }
    }
  }
  if (!inputDesignVec) {
    ierr = VecDestroy(&designVec);
    CHKERRQ(ierr);
  }
  if (!filter) {
    ierr = VecDestroy(&filteredVec);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
