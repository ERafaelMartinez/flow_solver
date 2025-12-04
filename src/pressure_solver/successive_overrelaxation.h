#pragma once

#include "../discretization/discretization.h"
#include "pressure_solver.h"
#include <iostream>
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

private:
    double omega_;  // Relaxation factor for SOR
};
