#include "simulation.h"
#include <cassert>
#include <cmath>
#ifndef DISABLE_OUTPUT_WRITERS
#include "../output_writer/output_writer_paraview_parallel.h"
#include "../output_writer/output_writer_text_parallel.h"
#endif
#include <algorithm>

// Constructor. Initializes the simulation, creating the staggered grid
// based on the simulation settings.
Simulation::Simulation(Settings *settings)
    : settings_(settings), time_step_(0.0), simulation_time_(0.0) {
  // Calculate cell size based on physical size and number of cells
  std::array<double, 2> cellSize = {
      settings->physicalSize[0] / settings->nCells[0],
      settings->physicalSize[1] / settings->nCells[1]};

#ifndef NDEBUG
  std::cout << "Initializing partitioning..." << std::endl;
#endif
  // Create and initialize domain partitioning
  partitioning_ = std::make_shared<Partitioning>();
  partitioning_->initialize(settings->nCells);

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] Initializing data exchanger..." << std::endl;
#endif
  // Initialize data exchanger
  dataExchanger_ = std::make_shared<DataExchanger>(partitioning_);

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] Initializing discretization..." << std::endl;
#endif
  // Initialize discretizations based on partitioning
  initDiscretization_(partitioning_->nCellsLocal(), cellSize);

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] Initializing pressure solver..." << std::endl;
#endif
  // Initialize pressure solver based on settings
  initPressureSolver_();

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] Initializing boundary manager..." << std::endl;
#endif
  // Initialize boundary manager
  initBoundaryManager_();

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] Initializing output writers..." << std::endl;
#endif
  // Initialize output writers
  initOutputWriters_();
}

// Create discretization based on settings
void Simulation::initDiscretization_(std::array<int, 2> nCells,
                                     std::array<double, 2> cellSize) {
  if (settings_->useDonorCell) {
#ifndef NDEBUG
    std::cout << "[" << partitioning_->ownRankNo() << "] Using donor cells!"
              << std::endl;
#endif
    discretization_ =
        std::make_shared<DonorCell>(nCells, cellSize, settings_->alpha);
  } else {
#ifndef NDEBUG
    std::cout << "[" << partitioning_->ownRankNo()
              << "] Using central differences!" << std::endl;
#endif
    discretization_ = std::make_shared<CentralDifferences>(nCells, cellSize);
  }

  // validate size of field variables:
  // each variable hast two ghost cells on each direction
  int nCellsExpected = (nCells[0] + 2) * (nCells[1] + 2);
  // assert(discretization_->u().data().size() == nCellsExpected);
  // assert(discretization_->v().data().size() == nCellsExpected);
}

// Create pressure solver based on settings
void Simulation::initPressureSolver_() {
  if (settings_->pressureSolver == "GaussSeidel") {
#ifndef NDEBUG
    std::cout << "[" << partitioning_->ownRankNo()
              << "] Using Gauss Seidel solver!" << std::endl;
#endif
    pressure_solver_ = std::make_shared<GaussSeidelPressureSolver>(
        discretization_, partitioning_, settings_->epsilon,
        settings_->maximumNumberOfIterations);
  } else if (settings_->pressureSolver == "SOR") {
#ifndef NDEBUG
    std::cout << "[" << partitioning_->ownRankNo() << "] Using SOR solver!"
              << std::endl;
#endif
    pressure_solver_ = std::make_shared<RedBlackSORPressureSolver>(
        discretization_, partitioning_, settings_->epsilon,
        settings_->maximumNumberOfIterations, settings_->omega);
  } else {
    throw std::invalid_argument("Unknown pressure solver type");
  }
}

// Initialize output writers
void Simulation::initOutputWriters_() {
// create writers
#ifndef DISABLE_OUTPUT_WRITERS
#ifndef NDEBUG
  writers_.push_back(std::make_unique<OutputWriterTextParallel>(discretization_,
                                                                partitioning_));
#endif
  writers_.push_back(std::make_unique<OutputWriterParaviewParallel>(
      discretization_, partitioning_));
#endif
}

// Initialize boundary manager
void Simulation::initBoundaryManager_() {
  boundaryManager_ = std::make_shared<BoundaryManager>(
      discretization_, partitioning_, settings_);
}

// Compute timestep based on the stability criteria
// derived from the grid size, reynolds number,
// and maximum velocities
double Simulation::computeNextTimeStepSize() {
  // difusion-based timestep restriction
  double dx = discretization_->cellSize()[0];
  double dy = discretization_->cellSize()[1];
  double diff_dt =
      0.5 * settings_->re * ((dx * dx * dy * dy) / (dx * dx + dy * dy));

  // convection-based timestep restriction

  // get the maximum value of the u and v fields (double-arrays)
  double u_max = discretization_->u().maxMagnitude();
  double v_max = discretization_->v().maxMagnitude();

  // exchange maximum (full domain) velocity with main rank
  std::array<double, 2> maxVelocity = {u_max, v_max};
  std::array<double, 2> globalMaxVelocity;
  globalMaxVelocity = dataExchanger_->getMaximumVelocity(maxVelocity);

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \t u_max, v_max = " << globalMaxVelocity[0] << ", "
            << globalMaxVelocity[1] << std::endl;
#endif

  double conv_dt_u = dx / std::abs(globalMaxVelocity[0]);
  double conv_dt_v = dy / std::abs(globalMaxVelocity[1]);
  double conv_dt = std::min(conv_dt_u, conv_dt_v);

  // take the smallest physics-induced dt and scale with safety factor
  double dt = std::min(diff_dt, conv_dt) * settings_->tau;

// limit the maximum timestep size using the settings
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \t Computed dt: " << std::min(dt, settings_->maximumDt)
            << std::endl;
#endif

  return std::min(dt, settings_->maximumDt);
}

// Compute the right-hand side for the pressure
// poisson equation
void Simulation::computeRHS() {
  // Implementation of RHS computation goes here
  double dx = discretization_->cellSize()[0];
  double dy = discretization_->cellSize()[1];
  double dt = time_step_;

  const FieldVariable &F = discretization_->f();
  const FieldVariable &G = discretization_->g();
  FieldVariable &rhs = discretization_->rhs();

  for (int j = discretization_->rhsJBegin(); j <= discretization_->rhsJEnd();
       j++) {
    for (int i = discretization_->rhsIBegin(); i <= discretization_->rhsIEnd();
         i++) {

      rhs.at(i, j) = (1. / dt) * (1 / dx * (F.at(i, j) - F.at(i - 1, j)) +
                                  1 / dy * (G.at(i, j) - G.at(i, j - 1)));
    }
  }
}

// Method to compute the intermediate velocities F, G
void Simulation::computeIntermediateVelocities() {
  // Implementation of intermediate velocity computation goes here
  double dt = time_step_;
  double Re = settings_->re;
  double gx = settings_->g[0];
  double gy = settings_->g[1];

  const FieldVariable &u = discretization_->u();
  const FieldVariable &v = discretization_->v();
  FieldVariable &F = discretization_->f();
  FieldVariable &G = discretization_->g();

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo() << "] \t\t Computing F..."
            << std::endl;
#endif

  for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); ++i) {
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd();
         ++j) {
      // Compute F(i,j)
      double d2udx2 = discretization_->computeD2uDx2(i, j);
      double d2udy2 = discretization_->computeD2uDy2(i, j);
      double du2dx = discretization_->computeDu2Dx(i, j);
      double duvdy = discretization_->computeDuvDy(i, j);
      double Aij = (1. / Re) * (d2udx2 + d2udy2) - du2dx - duvdy + gx;

      F.at(i, j) = u.at(i, j) + dt * Aij;
    }
  }

#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo() << "] \t\t Computing G..."
            << std::endl;
#endif
  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); ++i) {
    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd();
         ++j) {
      // Compute G(i,j)
      double d2vdx2 = discretization_->computeD2vDx2(i, j);
      double d2vdy2 = discretization_->computeD2vDy2(i, j);
      double dv2dy = discretization_->computeDv2Dy(i, j);
      double duvdx = discretization_->computeDuvDx(i, j);
      double Bij = (1. / Re) * (d2vdx2 + d2vdy2) - dv2dy - duvdx + gy;

      G.at(i, j) = v.at(i, j) + dt * Bij;
    }
  }
}

// Solve the pressure equation
void Simulation::solvePressureEquation() {
  pressure_solver_->solvePressureEquation();
}

// Recalculate the velocities based on the new pressure field
void Simulation::computeVelocities() {
  // Implementation of velocity computation goes here
  double dt = time_step_;

  FieldVariable &u = discretization_->u();
  FieldVariable &v = discretization_->v();
  const FieldVariable &F = discretization_->f();
  const FieldVariable &G = discretization_->g();

  for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); ++i) {
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd();
         ++j) {
      double dpdx = discretization_->computeDpDx(i, j);
      u.at(i, j) = F.at(i, j) - dt * dpdx;
    }
  }

  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); ++i) {
    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd();
         ++j) {
      double dpdy = discretization_->computeDpDy(i, j);
      v.at(i, j) = G.at(i, j) - dt * dpdy;
    }
  }
}

// Output the current state of the simulation
// using the OutputWritter class
void Simulation::outputSimulationState(double outputIndex) {
  for (std::vector< std::unique_ptr<OutputWriter> >::size_type i = 0;
       i < writers_.size(); i++) {
    writers_[i]->writeFile(outputIndex);
  }
}

// run simulation timestep
void Simulation::runTimestep() {
// 0. Apply/set boundary conditions for velocity field
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tSetting velocity boundaries..." << std::endl;
#endif
  boundaryManager_->setBoundaryConditionsVelocity();

// 1.1 Compute next time step size based on the values of
// the current velocity field and the stability criteria
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo() << "] \tComputing timestep..."
            << std::endl;
#endif
  time_step_ = computeNextTimeStepSize();
  simulation_time_ += time_step_;

// 2.1 Compute intermediate velocities F, G
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tComputing intermediate velocity field" << std::endl;
#endif
  computeIntermediateVelocities();
// 2.2 Enforce boundary conditions for F and G
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tSetting boundaries for F and G..." << std::endl;
#endif
  boundaryManager_->setBoundaryConditionsFG();
// 2.3 Exchange F and G with neighbors: not required

// 3.1 Compute RHS for pressure poisson equation
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo() << "] \tComputing rhs..."
            << std::endl;
#endif
  computeRHS();
// 3.2 Exchange RHS with neighbors: not required

// 4. Solve pressure equation
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tSolving pressure equation..." << std::endl;
#endif
  solvePressureEquation();

// 5. Compute velocities based on new pressure field
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tComputing velocities..." << std::endl;
#endif
  computeVelocities();

// 6. Exchange velocity field with neighbors
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tExchanging horizontal velocities..." << std::endl;
#endif
  dataExchanger_->exchange(discretization_->u());
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo()
            << "] \tExchanging vertical velocities..." << std::endl;
#endif
  dataExchanger_->exchange(discretization_->v());

// 7. Output current state of the simulation
#ifndef NDEBUG
  std::cout << "[" << partitioning_->ownRankNo() << "] \tWriting simulation at "
            << simulation_time_ << std::endl;
#endif
  // output simulation state if simulated time is next second
  if (std::floor(simulation_time_) !=
      std::floor(simulation_time_ + time_step_)) {
    outputSimulationState(simulation_time_);
  }
}

// run the full simulation
void Simulation::run() {
  int stepNumber = 0;
  while (simulation_time_ < settings_->endTime) {
#ifndef NDEBUG
    std::cout << "[" << partitioning_->ownRankNo() << "] Simulation step "
              << stepNumber << ":" << std::endl;
    std::cout << "[" << partitioning_->ownRankNo()
              << "] \t Current simulation time: " << simulation_time_
              << std::endl;
#endif
    runTimestep();
    stepNumber++;
  }
}
