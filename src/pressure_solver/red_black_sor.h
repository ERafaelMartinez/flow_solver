#pragma once

#include "../discretization/discretization.h"
#include "pressure_solver.h"
#include <iostream>
/* This file declares the SORPressureSolver class, which is a specific
   implementation of the PressureSolver interface. It uses the Successive
   Over-Relaxation (SOR) iterative method to solve the pressure Poisson equation
   in a numerical simulation framework.
*/
class RedBlackSORPressureSolver : public PressureSolver {
public:
  // Constructor
  RedBlackSORPressureSolver(std::shared_ptr<Discretization> discretization,
                            std::shared_ptr<Partitioning> partitioning,
                            double convergence_tol, int max_iterations,
                            double omega);

  /**
   * @brief Calculate the pressure for one iteration using SOR.
   *
   * This method performs the pressure calculation for a single iteration step
   * using the Successive Over-Relaxation (SOR) method.
   */
  void calcPressureIter() override;

private:
  enum GridColor { RED, BLACK };

  double omega_; // Relaxation factor for SOR
                 //

  void doColorIteration(GridColor color, const double coeff, const double idx2,
                        const double idy2);
};
