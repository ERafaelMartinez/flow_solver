#include "gauss_seidel.h"

/* This file implements the GaussSeidelPressureSolver class, which is a specific 
   implementation of the PressureSolver interface. It uses the Gauss-Seidel 
   iterative method to solve the pressure Poisson equation in a numerical 
   simulation framework.
*/

// Constructor
GaussSeidelPressureSolver::GaussSeidelPressureSolver(
    std::shared_ptr<Discretization> discretization,
    double convergence_tol,
    int max_iterations
) : PressureSolver(discretization, convergence_tol, max_iterations) {}

// Calculate pressure for one iteration
void GaussSeidelPressureSolver::calcPressureIter() {
    // Parameters required for the calculation
    auto [dx, dy] = discretization_->cellSize();
    auto [Nx, Ny] = discretization_->gridSize();
    double coeff = dx*dx * dy*dy / (2 * (dx*dx + dy*dy));

    res_ = 0.0; // Reset residual for this iteration
    for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); ++j)
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); ++i){
            discretization_->p().at(i, j) = coeff * (
                (discretization_->p().at(i-1, j) + discretization_->p().at(i+1, j)) / (dx*dx) +
                (discretization_->p().at(i, j-1) + discretization_->p().at(i, j+1)) / (dy*dy) -
                discretization_->rhs().at(i, j) 
            );

            // Add the squared difference to the residual times normalization factor
            double laplacian = (
                (discretization_->p().at(i+1, j) - 2*discretization_->p().at(i, j) + discretization_->p().at(i-1, j)) / (dx*dx) +
                (discretization_->p().at(i, j+1) - 2*discretization_->p().at(i, j) + discretization_->p().at(i, j-1)) / (dy*dy)
            );
            double diff = discretization_->rhs().at(i, j) - laplacian;
            res_ += diff * diff / (Nx * Ny);
        }

}
