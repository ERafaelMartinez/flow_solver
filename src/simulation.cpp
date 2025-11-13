#include "simulation.h"
#include "settings/settings.h"

// Constructor. Initializes the simulation, creating the staggered grid
// based on the simulation settings.
Simulation::Simulation(Settings* settings):
    settings_(settings),
    time_step_(0.0),
    simulation_time_(0.0)
{
    // Initialize staggered-grid/grid-discretization based on settings
    if (settings->useDonorCell) {
        discretization_ = new DonorCell(
            settings->nCells,
            {settings->physicalSize[0] / settings->nCells[0],
            settings->physicalSize[1] / settings->nCells[1]},
            settings->alpha);
        } 
    else {
    discretization_ = new CentralDifferences(
        settings->nCells,
        {settings->physicalSize[0] / settings->nCells[0],
        settings->physicalSize[1] / settings->nCells[1]
        });
    }


    // Initialize pressure solver based on settings
    if (settings->pressureSolver == "GaussSeidel") {
        pressure_solver_ = new GaussSeidelPressureSolver(
            discretization_,
            &settings->epsilon,
            &settings->maximumNumberOfIterations
        );
    }
    else if (settings->pressureSolver == "SOR") {
        throw std::invalid_argument("SOR solver not implemented yet");
    }
    else {
        throw std::invalid_argument("Unknown pressure solver type");
    }

}

// Destructor. Cleans up allocated resources.
Simulation::~Simulation() {
    delete grid_;
    delete discretization_;
    delete pressure_solver_;
}

// Method to apply/set boundary conditions for the velocity field
void Simulation::setBoundaryConditionsVelocity() {
    // Apply inhomogeneous Dirichlet boundary conditions
    // to the velocity u, and v fields based on settings
    // and apply them for F and G based on the Neumann BC for the
    // pressure posisson equation.

    // get reference to u, v, F, and G field variables
    FieldVariable u = discretization_->u();
    FieldVariable v = discretization_->v();
    FieldVariable f = discretization_->f();
    FieldVariable g = discretization_->g();

    // Apply top and bottom boundary conditions
    std::array<double, 2> dirichletBcBottom = settings_->dirichletBcBottom;
    std::array<double, 2> dirichletBcTop = settings_->dirichletBcTop;
    for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); ++i) {
        // Bottom boundary: u(i,0) = 2*u_bbc - u(i,1)
        // No point directly at the boundary for u, so we set the value
        // in the ghost cell so that the interpolation to the first
        // point in the grid yields the BC at the boundary.
        u(i, discretization_->uJBegin() - 1) = (
            2 * dirichletBcBottom[0] - u(i, discretization_->uJBegin())
        );
        // For v, as v is defined at the top/bottom edges of the cells
        // then we can directly set the value to the ghost cell.
        v(i, discretization_->vJBegin() - 1) = dirichletBcBottom[1];
        // G's boundary condition is derived/chosen from the Neumann BC for p:
        // g(i, 0) = v(i, 0) = dirichletBcBottom[1]
        g(i, discretization_->vJBegin() - 1) = dirichletBcBottom[1];

        // Top boundary: u(i,end) = 2*u_tbc - u(i,end-1)
        // Analogous as in the bottom case
        u(i, discretization_->uJEnd() + 1) = (
            2 * dirichletBcTop[0] - u(i, discretization_->uJEnd())
        );
        // For v there is a point directly at the boundary, so we set
        // it explicitly on the final valid point and propagate it to 
        // the ghost cell as well.
        v(i, discretization_->vJEnd()) = dirichletBcTop[1];
        v(i, discretization_->vJEnd() + 1) = dirichletBcTop[1];
        // G's boundary condition is derived/chosen from the Neumann BC for p:
        // g(i, jmax) = v(i, jmax) = dirichletBcTop[1]
        g(i, discretization_->vJEnd()) = dirichletBcTop[1];
        g(i, discretization_->vJEnd() + 1) = dirichletBcTop[1];
    }

    // Apply left and right boundary conditions
    std::array<double, 2> dirichletBcLeft = settings_->dirichletBcLeft;
    std::array<double, 2> dirichletBcRight = settings_->dirichletBcRight;
    // We start and end in the indices of the ghost cells to include the left-right
    // boundary conditions to the corners of the extended/ghost domain as well
    for (int j = discretization_->vJBegin() - 1; j <= discretization_->vJEnd() + 1; ++j) {
        // Left boundary:
        // The u value of the ghost cell falls directly at the boundary so we set it
        // directly
        u(discretization_->uIBegin() - 1, j) = dirichletBcLeft[0];
        // For v, the first domain point is not at the boundary, so we set the
        // value at the ghost cell such that the average yields the desired condition
        // inbetween: v(0, j) = 2*v_bbc - v(1,j)
        v(discretization_->vIBegin() - 1, j) = (
            2 * dirichletBcLeft[1] - v(discretization_->vIBegin(), j)
        );
        // F's boundary condition is derived/chosen from the Neumann BC for p:
        // f(0, j) = u(0, j) = dirichletBcLeft[0]
        f(discretization_->uIBegin() - 1, j) = dirichletBcLeft[0];

        // Right boundary:
        // Analogous to the left boundary for the u velocity,
        // the end point of the domain is directly at the boundary
        // so we set it directly and propagate to the ghost cell
        u(discretization_->uIEnd(), j) = dirichletBcRight[0];
        u(discretization_->uIEnd() + 1, j) = dirichletBcRight[0];
        // For v, we set the ghost cell value such that the average
        // with the last domain point yields the desired BC at the boundary
        // v(end, j) = 2*v_rbc - v(end-1,j
        v(discretization_->vIEnd() + 1, j) = (
            2 * dirichletBcRight[1] - v(discretization_->vIEnd(), j)
        );
        // F's boundary condition is derived/chosen from the Neumann BC for p:
        // f(imax, j) = u(imax, j) = dirichletBcRight[0];
        f(discretization_->uIEnd(), j) = dirichletBcRight[0];
        f(discretization_->uIEnd() + 1, j) = dirichletBcRight[0];
    }

}


// Method to apply/set boundary conditions for the pressure field
void Simulation::setBoundaryConditionsPressure() {
    // Derived from the Neumann boundary conditions for pressure:

    // Top/Bottom boundaries:
    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); ++i) {
        // Bottom boundary:   p(i,0) = p(i,1)
        discretization_->p(i, discretization_->pJBegin() - 1) = 
            discretization_->p(i, discretization_->pJBegin());
        // Top boundary:  p(i,jmax + 1) = p(i,jmax)
        discretization_->p(i, discretization_->pJEnd() + 1) = 
            discretization_->p(i, discretization_->pJEnd());
    }

    // Left/Right boundaries:
    for (int j = discretization_->pJBegin() - 1; j <= discretization_->pJEnd() + 1; ++j) {
        // Left boundary:   p(0,j) = p(1,j)
        discretization_->p(discretization_->pIBegin() - 1, j) = 
            discretization_->p(discretization_->pIBegin(), j);
        // Right boundary:  p(imax + 1,j) = p(imax,j)
        discretization_->p(discretization_->pIEnd() + 1, j) = 
            discretization_->p(discretization_->pIEnd(), j);
    }

}

// Compute timestep based on the stability criteria
// derived from the grid size, reynolds number,
// and maximum velocities
double Simulation::computeNextTimeStepSize() {
    // difusion-based timestep restriction
    double dx = discretization_->dx();
    double dy = discretization_->dy();
    double diff_dt = 0.5 * settings_->re * (dx * dx * dy * dy) / (dx * dx + dy * dy);

    // convection-based timestep restriction
    FieldVariable u_field = discretization_->u();
    FieldVariable v_field = discretization_->v();

    // get the maximum value of the u and v fields (double-arrays)
    double u_max = u_field.data().empty() ? 0.0 : 
        *std::max_element(
            discretization_->u().data().begin(),
            discretization_->u().data().end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); }
        );
    double v_max = discretization_->v().data().empty() ? 0.0 : 
        *std::max_element(
            discretization_->v().data().begin(),
            discretization_->v().data().end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); }
        );
    double conv_dt_u = dx / std::abs(u_max);
    double conv_dt_v = dy / std::abs(v_max);
    double conv_dt = std::min(conv_dt_u, conv_dt_v);
    
    return std::min(diff_dt, conv_dt) * settings_->tau;  // take the smallest and scale with safety factor
}

// Compute the right-hand side for the pressure 
// poisson equation
void Simulation::computeRHS() {
    // Implementation of RHS computation goes here
}

// Method to compute the intermediate velocities F, G
void Simulation::computeIntermediateVelocities() {
    // Implementation of intermediate velocity computation goes here
}

// Solve the pressure equation
void Simulation::solvePressureEquation() {
    pressure_solver_->solvePressureEquation();
}

// Recalculate the velocities based on the new pressure field
void Simulation::computeVelocities() {
    // Implementation of velocity computation goes here
}

// Output the current state of the simulation
// using the OutputWritter class
void Simulation::outputSimulationState(int outputIndex) {
    // Implementation of output writing goes here
}


// run simulation timestep
void Simulation::runTimestep(int stepNumber) {
    // 1. Apply boundary conditions for velocity, pressure,
    //    and intermediate velocity fields
    setBoundaryConditionsVelocity();
    setBoundaryConditionsPressure();

    // 2. Compute next time step size
    time_step_ = computeNextTimeStepSize();
    simulation_time_ += time_step_;

    // 3. Compute intermediate velocities F, G
    computeIntermediateVelocities();

    // 4. Compute RHS for pressure poisson equation
    computeRHS();

    // 5. Solve pressure equation
    solvePressureEquation();

    // 6. Compute velocities based on new pressure field
    computeVelocities();

    // 7. Output current state of the simulation
    outputSimulationState(stepNumber);
}

// run the full simulation
void Simulation::run() {
    int stepNumber = 0;
    while (simulation_time_ < settings_->endTime) {
        runTimestep(stepNumber);
        stepNumber++;
    }
}
