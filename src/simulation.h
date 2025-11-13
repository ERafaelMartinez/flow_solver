#pragma once

#include "discretization/discretization.h"
#include "discretization/central_differences.h"
#include "discretization/donor_cell.h"
#include "staggered_grid/staggered_grid.h"
#include "pressure_solver/pressure_solver.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/successive_overrelaxation.h"
#include "settings/settings.h"


class Simulation {
public:
    // Constructor
    Simulation(Settings* settings);

    // Destructor
    ~Simulation();

    // Method to apply/set boundary conditions for velocity field
    void setBoundaryConditionsVelocity();

    // Method to apply/set boundary conditions for pressure field
    void setBoundaryConditionsPressure();

    // Compute timestep based on the stability criteria
    // derived from the grid size, reynolds number,
    // and maximum velocities
    double computeNextTimeStepSize();

    // Compute the right-hand side for the pressure 
    // poisson equation
    void computeRHS();

    // Method to compute the intermediate velocities F, G
    void computeIntermediateVelocities();

    // Solve the pressure equation
    void solvePressureEquation();

    // Recalculate the velocities based on the new pressure field
    void computeVelocities();

    // Output the current state of the simulation
    // using the OutputWritter class
    void outputSimulationState();

    void runSimulationLoop();

private:
    Settings* settings_;
    StaggeredGrid* grid_;
    Discretization* discretization_;
    PressureSolver* pressure_solver_;

    double time_step_;
    double simulation_time_;

};
