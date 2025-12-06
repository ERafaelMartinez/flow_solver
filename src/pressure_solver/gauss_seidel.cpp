#include "gauss_seidel.h"

/* This file implements the GaussSeidelPressureSolver class, which is a specific
   implementation of the PressureSolver interface. It uses the Gauss-Seidel
   iterative method to solve the pressure Poisson equation in a numerical
   simulation framework.
*/

// Constructor
GaussSeidelPressureSolver::GaussSeidelPressureSolver(
    std::shared_ptr<Discretization> discretization,
    const Partitioning &partitioning, double convergence_tol,
    int max_iterations)
    : PressureSolver(discretization, partitioning, convergence_tol,
                     max_iterations) {}

// Calculate pressure for one iteration
void GaussSeidelPressureSolver::calcPressureIter() {
  // Parameters required for the calculation
  auto cellSize = discretization_->cellSize();
  double dx = cellSize[0];
  double dy = cellSize[1];
  auto gridSize = discretization_->gridSize();
  int Nx = gridSize[0];
  int Ny = gridSize[1];
  double coeff = dx * dx * dy * dy / (2 * (dx * dx + dy * dy));

  // Get references to pressure and rhs fields
  FieldVariable &p = discretization_->p();
  FieldVariable &rhs = discretization_->rhs();

  for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); ++i) {
    for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd();
         ++j) {
      p.at(i, j) =
          coeff * (
                      // p(i-1,j) (new) + p(i+1,j) (old)
                      (p.at(i - 1, j) + p.at(i + 1, j)) / (dx * dx) +
                      // p(i,j-1) (new) + p(i,j+1) (old)
                      (p.at(i, j - 1) + p.at(i, j + 1)) / (dy * dy) -
                      // rhs(i,j) (new; with velocities of previous timestep)
                      rhs.at(i, j));
    }
  }

  // Compute residual
  res_ = 0.; // Reset residual for this iteration
  for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); ++j)
    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd();
         ++i) {
      // Add the squared difference to the residual times normalization factor
      double laplacian =
          ((p.at(i + 1, j) - 2 * p.at(i, j) + p.at(i - 1, j)) / (dx * dx) +
           (p.at(i, j + 1) - 2 * p.at(i, j) + p.at(i, j - 1)) / (dy * dy));
      double diff = rhs.at(i, j) - laplacian;
      res_ += diff * diff / (Nx * Ny);
    }

  // at the end res_ = r**2 / (Nx*Ny), where r is the difference between rhs and
  // the new laplacian
}
