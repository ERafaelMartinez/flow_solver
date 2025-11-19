#include "pressure_solver.h"
#include <iostream>

/* This file implements the PressureSolver class, which is a specific
   implementation of the PressureSolver interface using the Gauss-Seidel
   iterative method to solve the pressure Poisson equation in a numerical
   simulation framework.
*/

// solve pressure equation until convergence
void PressureSolver::solvePressureEquation() {
  int iteration = 0;
  res_ = 1000.0; // Initialize with a large residual
  #ifndef NDEBUG
    std::cout << "\tStarting pressure solver..." << std::endl;
  #endif
  
  while (!solutionHasConverged() && iteration < max_iterations_) {
    calcPressureIter();

    #ifndef NDEBUG
      if (iteration % 100 == 0 || iteration == 0) {
          std::cout << "\t  Pressure solver iteration " << iteration << ", residual: " << res_ << std::endl;
      }
    #endif
    iteration++;
  }
  #ifndef NDEBUG
    if (solutionHasConverged()) {
        std::cout << "\tPressure solver converged in " << iteration << " iterations." << std::endl;
    } else {
        std::cout << "\tPressure solver did not converge in " << max_iterations_ << " iterations. Residual: " << res_ << std::endl;
    }
  #endif
}

// check for convergence
bool PressureSolver::solutionHasConverged() {
  // Assuming the res_ measure has been updated during the
  // last iteration
  return (res_ <= convergence_tol_ * convergence_tol_);
}
