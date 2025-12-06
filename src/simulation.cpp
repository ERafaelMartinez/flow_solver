#include "simulation.h"
#include "discretization/discretization.h"
#include <cassert>
#ifndef DISABLE_OUTPUT_WRITERS
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"
#endif
#include <algorithm>

// Constructor. Initializes the simulation, creating the staggered grid
// based on the simulation settings.
Simulation::Simulation(Settings *settings)
    : settings_(settings), time_step_(0.0), simulation_time_(0.0) {
  // Calculate cell size based on physical size and number of cells
  std::array<double, 2> cellSize = {
      settings->physicalSize[0] / settings->nCells[0],
      settings->physicalSize[1] / settings->nCells[1]
  };

  // Initialize domain partitioning
  partitioning_->initialize(settings->nCells);

  // Initialize data exchanger
  dataExchanger_ = std::make_shared<DataExchanger>(partitioning_);

  // Initialize discretizations based on partitioning
  initDiscretization_(partitioning_->nCellsLocal(), cellSize);

  // Initialize pressure solver based on settings
  initPressureSolver_();

  // Initialize output writers
  initOutputWriters_();
}

// Create discretization based on settings
void Simulation::initDiscretization_(
  std::array<int, 2> nCells, std::array<double, 2> cellSize) {
  if (settings_->useDonorCell) {
    #ifndef NDEBUG
        std::cout << "Using donor cells!" << std::endl;
    #endif
    discretization_ = std::make_shared<DonorCell>(nCells, cellSize, settings_->alpha);
  } else {
    #ifndef NDEBUG
        std::cout << "Using central differences!" << std::endl;
    #endif
    discretization_ = std::make_shared<CentralDifferences>(nCells, cellSize);
  }

  // validate size of field variables:
  // each variable hast two ghost cells on each direction
  int nCellsExpected = (nCells[0] + 2) * (nCells[1] + 2);
  assert(discretization_->u().data().size() == nCellsExpected);
  assert(discretization_->v().data().size() == nCellsExpected);
}

// Create pressure solver based on settings
void Simulation::initPressureSolver_() {
  if (settings_->pressureSolver == "GaussSeidel") {
    #ifndef NDEBUG
        std::cout << "Using Gauss Seidel solver!" << std::endl;
    #endif
    pressure_solver_ = std::make_shared<GaussSeidelPressureSolver>(
        discretization_, settings_->epsilon,
        settings_->maximumNumberOfIterations);
  } else if (settings_->pressureSolver == "SOR") {
    #ifndef NDEBUG
        std::cout << "Using SOR solver!" << std::endl;
    #endif
    pressure_solver_ = std::make_shared<SORPressureSolver>(
        discretization_, settings_->epsilon, settings_->maximumNumberOfIterations,
        settings_->omega);
  } else {
    throw std::invalid_argument("Unknown pressure solver type");
  }
}

// Initialize output writers
void Simulation::initOutputWriters_() {
  // create writers
  #ifndef DISABLE_OUTPUT_WRITERS
  #ifndef NDEBUG
    writers_.push_back(
      std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_)
    );
  #endif
    // BUG: Fix bug inside
    writers_.push_back(
      std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_)
    );
  #endif
}

// Destructor. Cleans up allocated resources.
Simulation::~Simulation() {
  // delete discretization_.get
  // delete pressure_solver_;
}

// Method to apply/set boundary conditions for the velocity field
void Simulation::setBoundaryConditionsVelocity() {
  // Apply inhomogeneous Dirichlet boundary conditions
  // to the velocity u, and v fields based on settings

  // get reference to u, v field variables
  FieldVariable &u = discretization_->u();
  FieldVariable &v = discretization_->v();

  // Apply top and bottom boundary conditions
  std::array<double, 2> dirichletBcBottom = settings_->dirichletBcBottom;
  std::array<double, 2> dirichletBcTop = settings_->dirichletBcTop;
  for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); ++i) {
    // No point directly at the boundary for u, so we set the value
    // in the ghost cell so that the interpolation to the first/end
    // point in the grid yields the BC at the boundary.

    // u(i,0) = 2*u_bbc - u(i, j_begin)
    u.at(i, discretization_->uJBegin() - 1) =
        (2 * dirichletBcBottom[0] - u.at(i, discretization_->uJBegin()));

    // Top boundary: u(i,end + 1) = 2*u_tbc - u(i,end)
    u.at(i, discretization_->uJEnd() + 1) =
        (2 * dirichletBcTop[0] - u.at(i, discretization_->uJEnd()));
  }

  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); ++i) {
    // For v, as v is defined at the top/bottom edges of the cells
    // then we can directly set the value to the ghost cell.
    v.at(i, discretization_->vJBegin() - 1) = dirichletBcBottom[1];

    // For v there is a point directly at the boundary, so we set
    // it explicitly on the final valid point and propagate it to
    // the ghost cell as well.
    v.at(i, discretization_->vJEnd()) = dirichletBcTop[1];
    // v.at(i, discretization_->vJEnd() + 1) = dirichletBcTop[1];
  }

  // Apply left and right boundary conditions
  std::array<double, 2> dirichletBcLeft = settings_->dirichletBcLeft;
  std::array<double, 2> dirichletBcRight = settings_->dirichletBcRight;

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions to the corners
  for (int j = discretization_->uJBegin() - 1;
       j <= discretization_->uJEnd() + 1; ++j) {
    // Left boundary:
    // The u value of the ghost cell falls directly at the boundary so we set it
    // directly
    // u(0, j) = u_lbc
    u.at(discretization_->uIBegin() - 1, j) = dirichletBcLeft[0];

    // Right boundary:
    // Analogous to the left boundary for the u velocity,
    // the end point of the domain is directly at the boundary
    // so we set it directly and propagate to the ghost cell
    // u(end, j) = u_rbc & u(end + 1, j) = u_rbc
    u.at(discretization_->uIEnd(), j) = dirichletBcRight[0];
    // u.at(discretization_->uIEnd() + 1, j) = dirichletBcRight[0];
  }

  for (int j = discretization_->vJBegin() - 1; j < discretization_->vJEnd() + 1;
       ++j) {
    // For v, the first domain point is not at the boundary, so we set the
    // value at the ghost cell such that the average yields the desired
    // condition inbetween: v(0, j) = 2*v_bbc - v(i_begin,j)
    v.at(discretization_->vIBegin() - 1, j) =
        (2 * dirichletBcLeft[1] - v.at(discretization_->vIBegin(), j));

    // For v, we set the ghost cell value such that the average
    // with the last domain point yields the desired BC at the boundary
    // v(end + 1, j) = 2*v_rbc - v(end,j)
    v.at(discretization_->vIEnd() + 1, j) =
        (2 * dirichletBcRight[1] - v.at(discretization_->vIEnd(), j));
  }
}

// Method to apply/set boundary conditions for F, and G
void Simulation::setBoundaryConditionsFG() {
  // Override Boundary values to F and G based to guarantee
  // Neumann BC for the pressure posisson equation.

  // get reference to f, g field variables
  FieldVariable &f = discretization_->f();
  FieldVariable &g = discretization_->g();
  const FieldVariable &v = discretization_->v();
  const FieldVariable &u = discretization_->u();

  // Apply top and bottom boundary conditions
  for (int i = discretization_->fIBegin(); i <= discretization_->fIEnd(); ++i) {
    // G's boundary condition is derived/chosen from the Neumann BC for p:
    // g(i, 0) = v(i, 0) = dirichletBcBottom[1]
    f.at(i, discretization_->fJBegin() - 1) =
        u.at(i, discretization_->fJBegin() - 1);

    // G's boundary condition is derived/chosen from the Neumann BC for p:
    // g(i, jmax) = v(i, jmax) = dirichletBcTop[1]
    f.at(i, discretization_->fJEnd() + 1) =
        u.at(i, discretization_->fJEnd() + 1);
    // g.at(i, discretization_->gJEnd() + 1) = v.at(i, discretization_->vJEnd()
    // + 1);
  }

  for (int i = discretization_->gIBegin(); i <= discretization_->gIEnd(); ++i) {
    g.at(i, discretization_->gJBegin() - 1) =
        v.at(i, discretization_->gJBegin() - 1);

    // F's boundary condition is derived/chosen from the Neumann BC for p:
    // f(i,jmax) = u(i,jmax) = dirichletBcTop[0];
    g.at(i, discretization_->gJEnd()) = v.at(i, discretization_->gJEnd());
  }

  // Apply left and right boundary conditions

  // We start and end in the indices of the ghost cells to include the
  // left-right boundary conditions to the corners of the extended/ghost domain
  // as well
  for (int j = discretization_->fJBegin() - 1;
       j <= discretization_->fJEnd() + 1; ++j) {

    // F's boundary condition is derived/chosen from the Neumann BC for p:
    // f(0, j) = u(0, j) = dirichletBcLeft[0]
    g.at(discretization_->gIBegin() - 1, j) =
        v.at(discretization_->gIBegin() - 1, j);

    // F's boundary condition is derived/chosen from the Neumann BC for p:
    // f(imax, j) = u(imax, j) = dirichletBcRight[0];
    // g.at(discretization_->gIEnd(), j) = v.at(discretization_->gIEnd(), j);
    g.at(discretization_->gIEnd() + 1, j) =
        v.at(discretization_->gIEnd() + 1, j);
  }

  for (int j = discretization_->fJBegin() - 1;
       j <= discretization_->fJEnd() + 1; ++j) {

    f.at(discretization_->fIBegin() - 1, j) =
        u.at(discretization_->fIBegin() - 1, j);

    f.at(discretization_->fIEnd(), j) = u.at(discretization_->fIEnd(), j);
  }
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

  #ifndef NDEBUG
    std::cout << "\t u_max, v_max = " << u_max << ", " << v_max << std::endl;
  #endif

  double conv_dt_u = dx / std::abs(u_max);
  double conv_dt_v = dy / std::abs(v_max);
  double conv_dt = std::min(conv_dt_u, conv_dt_v);

  // take the smallest physics-induced dt and scale with safety factor
  double dt = std::min(diff_dt, conv_dt) * settings_->tau;

  // limit the maximum timestep size using the settings
  #ifndef NDEBUG
    std::cout << "\t Computed dt: " << std::min(dt, settings_->maximumDt)
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
    std::cout << "\t\t Computing F..." << std::endl;
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
    std::cout << "\t\t Computing G..." << std::endl;
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

  for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); ++i) {
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd();
         ++j) {
      double dpdx = discretization_->computeDpDx(i, j);
      u.at(i, j) = F.at(i, j) - dt * dpdx;
    }
  }

  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); ++i) {
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd();
         ++j) {
      double dpdy = discretization_->computeDpDy(i, j);
      v.at(i, j) = G.at(i, j) - dt * dpdy;
    }
  }
}

// Output the current state of the simulation
// using the OutputWritter class
void Simulation::outputSimulationState(double outputIndex) {
  for (std::vector<std::unique_ptr<OutputWriter>>::size_type i = 0;
       i < writers_.size(); i++) {
    writers_[i]->writeFile(outputIndex);
  }
}

// run simulation timestep
void Simulation::runTimestep() {
  // 0. Apply/set boundary conditions for velocity field
  #ifndef NDEBUG
    std::cout << "\tSetting velocity boundaries..." << std::endl;
  #endif
  setBoundaryConditionsVelocity();

  // 1.1 Compute next time step size based on the values of
  // the current velocity field and the stability criteria
  #ifndef NDEBUG
    std::cout << "\tComputing timestep..." << std::endl;
  #endif
  time_step_ = computeNextTimeStepSize();
  // 1.2 obtain maximum time step size from main rank
  #ifndef NDEBUG
    std::cout << "\tExchanging timestep..." << std::endl;
  #endif
  time_step_ = dataExchanger_->getMaximumTimeStepSize(time_step_);
  simulation_time_ += time_step_;

  // 3. Compute intermediate velocities F, G
  #ifndef NDEBUG
    std::cout << "\tComputing intermediate velocity field" << std::endl;
  #endif
  computeIntermediateVelocities();

  // 2. Enforce boundary conditions for F and G
  #ifndef NDEBUG
    std::cout << "\tSetting boundaries for F and G..." << std::endl;
  #endif
  setBoundaryConditionsFG();

  // 4. Compute RHS for pressure poisson equation
  #ifndef NDEBUG
    std::cout << "\tComputing rhs..." << std::endl;
  #endif
  computeRHS();

  // 5. Solve pressure equation
  #ifndef NDEBUG
    std::cout << "\tSolving pressure equation..." << std::endl;
  #endif
  solvePressureEquation();

  // 6. Compute velocities based on new pressure field
  #ifndef NDEBUG
    std::cout << "\tComputing velocities..." << std::endl;
  #endif
  computeVelocities();

  // 7. Output current state of the simulation
  #ifndef NDEBUG
    std::cout << "\tWriting simulation at " << simulation_time_ << std::endl;
  #endif
  outputSimulationState(simulation_time_);
}

// run the full simulation
void Simulation::run() {
  int stepNumber = 0;
  while (simulation_time_ < settings_->endTime) {
    #ifndef NDEBUG
      std::cout << "Simulation step " << stepNumber << ":" << std::endl;
      std::cout << "\t Current simulation time: " << simulation_time_
                << std::endl;
    #endif
    runTimestep();
    stepNumber++;
  }
}
