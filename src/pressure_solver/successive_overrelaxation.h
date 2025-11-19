#pragma once

#include "discretization/discretization.h"
#include "pressure_solver.h"

/* This file declares the SORPressureSolver class, which is a specific implementation of the 
   PressureSolver interface. It uses the Successive Over-Relaxation (SOR) iterative method 
   to solve the pressure Poisson equation in a numerical simulation framework.
*/
class SORPressureSolver : public PressureSolver {
public:
    // Constructor
    SORPressureSolver(std::shared_ptr<Discretization> discretization,
        double convergence_tol,
        int max_iterations,
        double omega
    );

    /**
     * @brief Calculate the pressure for one iteration using SOR.
     * 
     * This method performs the pressure calculation for a single iteration step
     * using the Successive Over-Relaxation (SOR) method.
     */
    void calcPressureIter() override;

    /**
     * @brief Solve the pressure equation using SOR.
     * 
     * This method iteratively solves the pressure equation using the SOR method
     * until the solution converges or the maximum number of iterations is reached.
     */
    void solvePressureEquation();

    /**
     * @brief Check if the solution has converged in SOR.
     * 
     * This method evaluates whether the solution meets the convergence criteria
     * specific to the SOR method.
     * @return True if the solution has converged, false otherwise.
     */
    bool solutionHasConverged();

private:
    double omega_;  // Relaxation factor for SOR
};
