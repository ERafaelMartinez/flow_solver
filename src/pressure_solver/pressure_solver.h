#pragma once

//In this file, we define the PressureSolver class interface
#include "../discretization/discretization.h"
#include <memory>


/* This class defines the interface for pressure solvers in the numerical simulation framework.
   It optionally receives as input a discretization scheme derived from SystemDiscretization
   and provides methods to set the discretization and solve the pressure equation accordingly.
*/
class PressureSolver{
public:
    // Constructor
    PressureSolver(
        std::shared_ptr<Discretization> discretization,
        double convergence_tol,
        int max_iterations
    ) : 
    discretization_(discretization),
    convergence_tol_(convergence_tol), 
    max_iterations_(max_iterations) {}

   virtual ~PressureSolver() {}

    /**
     * @brief Calculate the pressure for one iteration.
     * 
     * This is a pure virtual method that must be implemented by derived classes.
     * It performs the pressure calculation for a single iteration step.
     */
    virtual void calcPressureIter() = 0;

    /**
     * @brief Perform multiple iterations until convergence.
     * 
     * This method iteratively solves the pressure equation until the solution
     * converges or the maximum number of iterations is reached.
     */
    void solvePressureEquation();

    void calcRes();

    /**
     * @brief Check if the solution has converged.
     * 
     * This method evaluates whether the solution meets the convergence criteria.
     * @return True if the solution has converged, false otherwise.
     */
    bool solutionHasConverged();

    /**
     * @brief Set the boundary conditions for the pressure field.
     * 
     * This method applies the appropriate homogeneous Neumann boundary conditions
     * to the pressure field in the discretization.
     */
    void setBoundaryConditions();

protected:
    std::shared_ptr<Discretization> discretization_;
    double res_;  // Residual measure for convergence check, has to be calculated/updated in calcPressureIter
    double convergence_tol_;
    int max_iterations_;
};
