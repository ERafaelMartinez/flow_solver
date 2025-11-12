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

    // Method to calculate the pressure for one iteration
    virtual void calcPressureIter() = 0;

    // Method to perform multiple iterations until convergence
    void solvePressureEquation();

    // Method to check for convergence
    bool solutionHasConverged();

protected:
    Discretization* discretization_;
    double* res_;  // Residual measure for convergence check, has to be calculated/updated in calcPressureIter
    double* convergence_tol_;
    int* max_iterations_;
};
