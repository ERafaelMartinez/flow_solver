#pragma once

#include "discretization/discretization.h"
#include "pressure_solver.h"

/* This file declares the GaussSeidelPressureSolver class, which is a specific 
   implementation of the PressureSolver interface. It uses the Gauss-Seidel 
   iterative method to solve the pressure Poisson equation in a numerical 
   simulation framework.
*/

class GaussSeidelPressureSolver : public PressureSolver {
public:
    // Constructor
    GaussSeidelPressureSolver(
        std::shared_ptr<Discretization> discretization,
        double convergence_tol,
        int max_iterations
    );

    // Override methods from PressureSolver
    /**
     * @brief Calculate the pressure for one iteration using Gauss-Seidel.
     * 
     * This method performs the pressure calculation for a single iteration step
     * using the Gauss-Seidel method.
     */
    void calcPressureIter() override;

    /**
     * @brief Solve the pressure equation using Gauss-Seidel.
     * 
     * This method iteratively solves the pressure equation using the Gauss-Seidel method
     * until the solution converges or the maximum number of iterations is reached.
     */
    void solvePressureEquation();

    /**
     * @brief Check if the solution has converged in Gauss-Seidel.
     * 
     * This method evaluates whether the solution meets the convergence criteria
     * specific to the Gauss-Seidel method.
     * @return True if the solution has converged, false otherwise.
     */
    bool solutionHasConverged();
private:
    // Private members specific to Gauss-Seidel implementation
};
