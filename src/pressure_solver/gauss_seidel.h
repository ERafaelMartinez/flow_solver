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
        Discretization* discretization,
        double* convergence_tol_,
        int* max_iterations_
    );

    // Override methods from PressureSolver
    void calcPressureIter() override;
    void solvePressureEquation() override;
    bool solutionHasConverged() override;
private:
    // Private members specific to Gauss-Seidel implementation
};  
