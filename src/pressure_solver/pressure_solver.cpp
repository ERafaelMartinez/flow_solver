#include "pressure_solver.h"
#include <iostream>
#include <ostream>
#include <cmath>

/* This file implements the PressureSolver class, which is a specific
   implementation of the PressureSolver interface using the Gauss-Seidel
   iterative method to solve the pressure Poisson equation in a numerical
   simulation framework.
*/

// solve pressure equation until convergence
void PressureSolver::solvePressureEquation() {
  int iteration = 0;
  res_ = 1000.0; // Initialize with a large residual
#ifndef NDEBUG
  std::cout << "\tStarting pressure solver..." << std::endl;
#endif

  setBoundaryConditions();

  while (!solutionHasConverged() && iteration < max_iterations_) {
    calcPressureIter();
    dataExchanger_->exchange(discretization_->p());
    setBoundaryConditions();
    calcRes();
    iteration++;
  }
#ifndef NDEBUG
  if (solutionHasConverged()) {
    std::cout << "\tPressure solver converged in " << iteration
              << " iterations." << std::endl;
  } else {
    std::cout << "\tPressure solver did not converge in " << max_iterations_
              << " iterations. Residual: " << res_ << std::endl;
  }
#endif
}

void PressureSolver::calcRes() {
  double localResidual = 0.0;

  auto cellSize = discretization_->cellSize();
  const double dx = cellSize[0];
  const double dy = cellSize[1];
  const double idx2 = 1.0 / (dx * dx);
  const double idy2 = 1.0 / (dy * dy);

  const int Nx = discretization_->gridSize()[0];
  const int Ny = discretization_->gridSize()[1];

  FieldVariable &p = discretization_->p();
  const FieldVariable &rhs = discretization_->rhs();

  for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
    for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd();
         j++) {
      double laplacian =
          ((p.at(i + 1, j) - 2 * p.at(i, j) + p.at(i - 1, j)) * idx2 +
           (p.at(i, j + 1) - 2 * p.at(i, j) + p.at(i, j - 1)) * idy2);
      double diff = std::fabs(laplacian - rhs.at(i, j));
      localResidual += (diff * diff);
    }
  }

  localResidual /= (Nx * Ny);

  // obtain total residual from all ranks
  res_ = dataExchanger_->getResidual(localResidual);
}

// check for convergence
bool PressureSolver::solutionHasConverged() {
  // Assuming the res_ measure has been updated during the
  // last iteration
  return (res_ <= convergence_tol_ * convergence_tol_);
}

// Method to apply/set boundary conditions for the pressure field
void PressureSolver::setBoundaryConditions() {
  if (partitioning_->ownPartitionContainsBottomBoundary()) {
    setBoundaryConditionsBottom();
  }
  if (partitioning_->ownPartitionContainsTopBoundary()) {
    setBoundaryConditionsTop();
  }
  if (partitioning_->ownPartitionContainsLeftBoundary()) {
    setBoundaryConditionsLeft();
  }
  if (partitioning_->ownPartitionContainsRightBoundary()) {
    setBoundaryConditionsRight();
  }
}

// Apply bottom boundary conditions
void PressureSolver::setBoundaryConditionsBottom() {
  FieldVariable &pressure_field = discretization_->p();

  // Bottom boundary: p(i,0) = p(i,1)
  for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
    pressure_field.at(i, discretization_->pJBegin() - 1) =
        pressure_field.at(i, discretization_->pJBegin());
  }
}

// Apply top boundary conditions
void PressureSolver::setBoundaryConditionsTop() {
  FieldVariable &pressure_field = discretization_->p();

  // Top boundary: p(i,jmax + 1) = p(i,jmax)
  for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
    pressure_field.at(i, discretization_->pJEnd() + 1) =
        pressure_field.at(i, discretization_->pJEnd());
  }
}

// Apply left boundary conditions
void PressureSolver::setBoundaryConditionsLeft() {
  FieldVariable &pressure_field = discretization_->p();

  // Left boundary: p(0,j) = p(1,j)
  for (int j = discretization_->pJBegin() - 1; j <= discretization_->pJEnd() + 1; j++) {
    pressure_field.at(discretization_->pIBegin() - 1, j) =
        pressure_field.at(discretization_->pIBegin(), j);
  }
}

// Apply right boundary conditions
void PressureSolver::setBoundaryConditionsRight() {
  FieldVariable &pressure_field = discretization_->p();

  // Right boundary: p(imax + 1,j) = p(imax,j)
  for (int j = discretization_->pJBegin() - 1; j <= discretization_->pJEnd() + 1; j++) {
    pressure_field.at(discretization_->pIEnd() + 1, j) =
        pressure_field.at(discretization_->pIEnd(), j);
  }
}
