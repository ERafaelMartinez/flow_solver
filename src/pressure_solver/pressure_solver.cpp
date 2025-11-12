#include "pressure_solver.h"

/* This file implements the PressureSolver class, which is a specific
   implementation of the PressureSolver interface using the Gauss-Seidel iterative
   method to solve the pressure Poisson equation in a numerical simulation framework.
*/

// Constructor
PressureSolver::PressureSolver(
      Discretization* discretization,
      double* convergence_tol_,
      int* max_iterations_
   ) :
   discretization_(discretization),
   convergence_tol_(convergence_tol_),
   max_iterations_(max_iterations_) {}

// Destructor
PressureSolver::~PressureSolver() {}

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
