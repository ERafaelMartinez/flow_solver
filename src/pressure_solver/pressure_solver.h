//In this file, we define the PressureSolver class interface
#include "discretization/discretization.h"
#include "staggered_grid/field_variable.h"
#include "settings/settings.h"


/* This class defines the interface for pressure solvers in the numerical simulation framework.
   It optionally receives as input a discretization scheme derived from SystemDiscretization
   and provides methods to set the discretization and solve the pressure equation accordingly.
*/
class PressureSolver{
public:
    // Constructor
    PressureSolver(
        Discretization* discretization,
        double* convergence_tol_,
        int* max_iterations_
    ) : 
    discretization_(discretization),
    convergence_tol_(convergence_tol_), 
    max_iterations_(max_iterations_) {}

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

    /**
     * @brief Check if the solution has converged.
     * 
     * This method evaluates whether the solution meets the convergence criteria.
     * @return True if the solution has converged, false otherwise.
     */
    bool solutionHasConverged();

protected:
    Discretization* discretization_;
    double* res_;  // Residual measure for convergence check, has to be calculated/updated in calcPressureIter
    double* convergence_tol_;
    int* max_iterations_;
};
