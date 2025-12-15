#include "gauss_seidel.h"

/* This file implements the GaussSeidelPressureSolver class, which is a specific
   implementation of the PressureSolver interface. It uses the Gauss-Seidel
   iterative method to solve the pressure Poisson equation in a numerical
   simulation framework.
*/

// Constructor
GaussSeidelPressureSolver::GaussSeidelPressureSolver(
    std::shared_ptr<Discretization> discretization,
    std::shared_ptr<Partitioning> partitioning, double convergence_tol,
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
  const double idx2 = 1.0 / (dx * dx);
  const double idy2 = 1.0 / (dy * dy);
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
                      (p.at(i - 1, j) + p.at(i + 1, j)) * idx2 +
                      // p(i,j-1) (new) + p(i,j+1) (old)
                      (p.at(i, j - 1) + p.at(i, j + 1)) * idy2 -
                      // rhs(i,j) (new; with velocities of previous timestep)
                      rhs.at(i, j));
    }
  }

}
