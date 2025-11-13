#include "pressure_solver.h"

/* This file implements the PressureSolver class, which is a specific
   implementation of the PressureSolver interface using the Gauss-Seidel iterative
   method to solve the pressure Poisson equation in a numerical simulation framework.
*/

// solve pressure equation until convergence
void PressureSolver::solvePressureEquation() {
    int iteration = 0;
    while (!solutionHasConverged() && iteration < *max_iterations_) {
        calcPressureIter();
        iteration++;
    }
}

// check for convergence
bool PressureSolver::solutionHasConverged() {
    // Assuming the res_ measure has been updated during the 
    // last iteration
   return (*res_ < *convergence_tol_*(*convergence_tol_));

}
