#pragma once

#include "discretization/central_differences.h"
#include "discretization/discretization.h"
#include "discretization/donor_cell.h"
#include "output_writer/output_writer.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/pressure_solver.h"
#include "pressure_solver/successive_overrelaxation.h"
#include "settings/settings.h"
#include "staggered_grid/staggered_grid.h"
#include <memory>
#include <vector>

class Simulation {
public:
  // Constructor
  Simulation(Settings *settings);

  // Destructor
  ~Simulation();

  // Method to apply/set boundary conditions for velocity field
  void setBoundaryConditionsVelocity();

  // Compute timestep based on the stability criteria
  // derived from the grid size, reynolds number,
  // and maximum velocities
  double computeNextTimeStepSize();

  // Compute the right-hand side for the pressure
  // poisson equation
  void computeRHS();

  // Method to compute the intermediate velocities F, G
  void computeIntermediateVelocities();

  // Method to apply/set boundary conditions for F and G
  void setBoundaryConditionsFG();

  // Solve the pressure equation
  void solvePressureEquation();

  // Recalculate the velocities based on the new pressure field
  void computeVelocities();

  // Output the current state of the simulation
  // using the OutputWritter class
  void outputSimulationState(double outputIndex);

  void runTimestep();

  void run();

private:
  Settings *settings_;
  std::shared_ptr<Partitioning> partitioning_;
  std::shared_ptr<Discretization> discretization_;
  std::shared_ptr<PressureSolver> pressure_solver_;
  std::vector< std::unique_ptr<OutputWriter> > writers_;

  void initPartitioning_();

  void initDiscretization_(std::array<int, 2> nCells, std::array<double, 2> cellSize);

  void initPressureSolver_();

  void initOutputWriters_();

  double time_step_;
  double simulation_time_;
};
