#include "successive_overrelaxation.h"

/* This file implements the SORPressureSolver class, which is a specific 
   implementation of the PressureSolver interface. It uses the SOR 
   iterative method to solve the pressure Poisson equation in a numerical 
   simulation framework.
*/

// Constructor
SORPressureSolver::SORPressureSolver(
    std::shared_ptr<Discretization> discretization,
    double convergence_tol,
    int max_iterations,
    double omega
) : PressureSolver(discretization, convergence_tol, max_iterations),
    omega_(omega)
{}

// Calculate pressure for one iteration
void SORPressureSolver::calcPressureIter() {
    // Parameters required for the calculation
    auto [dx, dy] = discretization_->cellSize();
    auto [Nx, Ny] = discretization_->gridSize();
    
    double coeff = dx*dx * dy*dy / (2 * (dx*dx + dy*dy));
    res_ = 0.0;

    for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); ++j)
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); ++i) {

            // Gauss-Seidel new pressure (unrelaxed)
            double p_GS = coeff * (
                (discretization_->p().at(i-1, j) + discretization_->p().at(i+1, j)) / (dx*dx) +
                (discretization_->p().at(i, j-1) + discretization_->p().at(i, j+1)) / (dy*dy) -
                discretization_->rhs().at(i, j)
            );

            // Apply SOR relaxation
            double p_old = discretization_->p().at(i, j);
            double p_relaxed = (1.0 - omega_) * p_old + omega_ * p_GS;

            discretization_->p().at(i, j) = p_relaxed;

            // Residual computation
            double laplacian =
                (discretization_->p().at(i+1, j) - 2*p_relaxed + discretization_->p().at(i-1, j)) / (dx*dx) +
                (discretization_->p().at(i, j+1) - 2*p_relaxed + discretization_->p().at(i, j-1)) / (dy*dy);

            double diff = discretization_->rhs().at(i, j) - laplacian;
            res_ += diff * diff / (Nx * Ny);
        }
}