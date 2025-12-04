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
    auto cellSize = discretization_->cellSize();
    const double dx = cellSize[0];
    const double dy = cellSize[1];
    const double idx2 = 1.0 / (dx * dx);
    const double idy2 = 1.0 / (dy * dy);
    
    // Precompute coefficient
    const double coeff = 1.0 / (2.0 * (idx2 + idy2));
    
    FieldVariable& p = discretization_->p();
    const FieldVariable& rhs = discretization_->rhs();
    
    for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
            // Compute the Gauss-Seidel target value
            const double p_GS = coeff * (
                (p.at(i-1, j) + p.at(i+1, j)) * idx2 +
                (p.at(i, j-1) + p.at(i, j+1)) * idy2 -
                rhs.at(i, j)
            );

            const double p_old = p.at(i, j);
            const double p_new = (1.0 - omega_) * p_old + omega_ * p_GS;
            p.at(i, j) = p_new;
        }
    }
}