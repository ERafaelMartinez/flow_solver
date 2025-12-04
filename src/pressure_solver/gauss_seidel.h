#pragma once

#include "../discretization/discretization.h"
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
private:
    // Private members specific to Gauss-Seidel implementation
};
