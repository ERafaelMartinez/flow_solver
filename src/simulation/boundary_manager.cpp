#include "boundary_manager.h"
#include "../storage/field_variable.h"

BoundaryManager::BoundaryManager(std::shared_ptr<Discretization> discretization,
                                 std::shared_ptr<Partitioning> partitioning,
                                 Settings *settings)
    : discretization_(discretization), partitioning_(partitioning),
      settings_(settings) {}

// =============================================================================
// Public Methods
// =============================================================================

void BoundaryManager::setBoundaryConditionsVelocity() {
  // Apply inhomogeneous Dirichlet boundary conditions
  // to the velocity u and v fields based on settings

  // The boundary conditions are only applied to the bounds of
  // the subdomain which correspond to the domain boundaries

  // Apply each boundary side separately for parallel control
  if (partitioning_->ownPartitionContainsTopBoundary()) {
    applyTopBoundaryU();
    applyTopBoundaryV();
  }

  if (partitioning_->ownPartitionContainsBottomBoundary()) {
    applyBottomBoundaryU();
    applyBottomBoundaryV();
  }

  if (partitioning_->ownPartitionContainsLeftBoundary()) {
    applyLeftBoundaryU();
    applyLeftBoundaryV();
  }

  if (partitioning_->ownPartitionContainsRightBoundary()) {
    applyRightBoundaryU();
    applyRightBoundaryV();
  }
}

void BoundaryManager::setBoundaryConditionsFG() {
  // Override boundary values of F and G to guarantee
  // Neumann BC for the pressure Poisson equation.

  // The boundary conditions are only applied to the bounds of
  // the subdomain which correspond to the domain boundaries

  // Apply each boundary side separately for parallel control
    applyTopBoundaryF();
    applyTopBoundaryG();

    applyBottomBoundaryF();
    applyBottomBoundaryG();

    applyLeftBoundaryF();
    applyLeftBoundaryG();

    applyRightBoundaryF();
    applyRightBoundaryG();
}

// =============================================================================
// Private Methods - Velocity Boundary Conditions (u, v)
// =============================================================================

void BoundaryManager::applyTopBoundaryU() {
  FieldVariable &u = discretization_->u();
  const std::array<double, 2> &dirichletBcTop = settings_->dirichletBcTop;

  for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); ++i) {
    // No point directly at the top boundary for u, so we set the value
    // in the ghost cell so that the interpolation to the end
    // point in the grid yields the BC at the boundary.

    // Top boundary: u(i,end + 1) = 2*u_tbc - u(i,end)
    u.at(i, discretization_->uJEnd() + 1) =
        (2 * dirichletBcTop[0] - u.at(i, discretization_->uJEnd()));
  }
}

void BoundaryManager::applyBottomBoundaryU() {
  FieldVariable &u = discretization_->u();
  const std::array<double, 2> &dirichletBcBottom = settings_->dirichletBcBottom;

  for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); ++i) {
    // No point directly at the boundary for u, so we set the value
    // in the ghost cell so that the interpolation to the first
    // point in the grid yields the BC at the boundary.

    // Bottom boundary: u(i,0) = 2*u_bbc - u(i, j_begin)
    u.at(i, discretization_->uJBegin() - 1) =
        (2 * dirichletBcBottom[0] - u.at(i, discretization_->uJBegin()));
  }
}

void BoundaryManager::applyTopBoundaryV() {
  FieldVariable &v = discretization_->v();
  const std::array<double, 2> &dirichletBcTop = settings_->dirichletBcTop;

  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); ++i) {
    // For v there is a point directly at the boundary, so we set
    // it explicitly on the final valid point and propagate it to
    // the ghost cell as well.
    v.at(i, discretization_->vJEnd()) = dirichletBcTop[1];
  }
}

void BoundaryManager::applyBottomBoundaryV() {
  FieldVariable &v = discretization_->v();
  const std::array<double, 2> &dirichletBcBottom = settings_->dirichletBcBottom;

  for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); ++i) {
    // For v, as v is defined at the top/bottom edges of the cells
    // then we can directly set the value to the ghost cell.
    v.at(i, discretization_->vJBegin() - 1) = dirichletBcBottom[1];
  }
}

void BoundaryManager::applyLeftBoundaryU() {
  FieldVariable &u = discretization_->u();
  const std::array<double, 2> &dirichletBcLeft = settings_->dirichletBcLeft;

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions over the corners
  for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd() + 1; ++j) {
    // Left boundary:
    // The u value of the ghost cell falls directly at the boundary so we set it
    // directly
    // u(0, j) = u_lbc
    u.at(discretization_->uIBegin() - 1, j) = dirichletBcLeft[0];
  }
}

void BoundaryManager::applyRightBoundaryU() {
  FieldVariable &u = discretization_->u();
  const std::array<double, 2> &dirichletBcRight = settings_->dirichletBcRight;

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions over the corners
  for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd() + 1; ++j) {
    // Right boundary:
    // Analogous to the left boundary for the u velocity,
    // the end point of the domain is directly at the boundary
    // so we set it directly and propagate to the ghost cell
    // u(end, j) = u_rbc & u(end + 1, j) = u_rbc
    u.at(discretization_->uIEnd(), j) = dirichletBcRight[0];
  }
}

void BoundaryManager::applyLeftBoundaryV() {
  FieldVariable &v = discretization_->v();
  const std::array<double, 2> &dirichletBcLeft = settings_->dirichletBcLeft;

  for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); ++j) {
    // For v, the first domain point is not at the boundary, so we set the
    // value at the ghost cell such that the average yields the desired
    // condition inbetween: v(0, j) = 2*v_bbc - v(i_begin,j)
    v.at(discretization_->vIBegin() - 1, j) =
        (2 * dirichletBcLeft[1] - v.at(discretization_->vIBegin(), j));
  }
}

void BoundaryManager::applyRightBoundaryV() {
  FieldVariable &v = discretization_->v();
  const std::array<double, 2> &dirichletBcRight = settings_->dirichletBcRight;

  for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); ++j) { 
    // For v, we set the ghost cell value such that the average
    // with the last domain point yields the desired BC at the boundary
    // v(end + 1, j) = 2*v_rbc - v(end,j)
    v.at(discretization_->vIEnd() + 1, j) =
        (2 * dirichletBcRight[1] - v.at(discretization_->vIEnd(), j));
  }
}

// =============================================================================
// Private Methods - Intermediate Velocity Boundary Conditions (F, G)
// =============================================================================

void BoundaryManager::applyTopBoundaryF() {
  FieldVariable &f = discretization_->f();
  const FieldVariable &u = discretization_->u();

  for (int i = discretization_->fIBegin(); i <= discretization_->fIEnd(); ++i) {
    // F's boundary condition is derived from the Neumann BC for p:
    // f(i,jmax + 1) = u(i,jmax + 1)
    f.at(i, discretization_->fJEnd() + 1) =
        u.at(i, discretization_->fJEnd() + 1);
  }
}

void BoundaryManager::applyBottomBoundaryF() {
  FieldVariable &f = discretization_->f();
  const FieldVariable &u = discretization_->u();

  for (int i = discretization_->fIBegin(); i <= discretization_->fIEnd(); ++i) {
    // F's boundary condition is derived from the Neumann BC for p:
    // f(i, 0) = u(i, 0)
    f.at(i, discretization_->fJBegin() - 1) =
        u.at(i, discretization_->fJBegin() - 1);
  }
}

void BoundaryManager::applyTopBoundaryG() {
  FieldVariable &g = discretization_->g();
  const FieldVariable &v = discretization_->v();

  for (int i = discretization_->gIBegin(); i <= discretization_->gIEnd(); ++i) {
    // G's boundary condition is derived from the Neumann BC for p:
    // g(i, jmax) = v(i, jmax)
    g.at(i, discretization_->gJEnd()) = v.at(i, discretization_->gJEnd());
  }
}

void BoundaryManager::applyBottomBoundaryG() {
  FieldVariable &g = discretization_->g();
  const FieldVariable &v = discretization_->v();

  for (int i = discretization_->gIBegin(); i <= discretization_->gIEnd(); ++i) {
    // G's boundary condition is derived from the Neumann BC for p:
    // g(i, 0) = v(i, 0)
    g.at(i, discretization_->gJBegin() - 1) =
        v.at(i, discretization_->gJBegin() - 1);
  }
}

void BoundaryManager::applyLeftBoundaryF() {
  FieldVariable &f = discretization_->f();
  const FieldVariable &u = discretization_->u();

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions over the corners
  for (int j = discretization_->fJBegin(); j <= discretization_->fJEnd(); ++j) {
    // F's boundary condition is derived from the Neumann BC for p:
    // f(0, j) = u(0, j)
    f.at(discretization_->fIBegin() - 1, j) =
        u.at(discretization_->fIBegin() - 1, j);
  }
}

void BoundaryManager::applyRightBoundaryF() {
  FieldVariable &f = discretization_->f();
  const FieldVariable &u = discretization_->u();

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions over the corners
  for (int j = discretization_->fJBegin(); j <= discretization_->fJEnd(); ++j) {
    // F's boundary condition is derived from the Neumann BC for p:
    // f(imax, j) = u(imax, j)
    f.at(discretization_->fIEnd(), j) = u.at(discretization_->fIEnd(), j);
  }
}

void BoundaryManager::applyLeftBoundaryG() {
  FieldVariable &g = discretization_->g();
  const FieldVariable &v = discretization_->v();

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions over the corners
  for (int j = discretization_->gJBegin(); j <= discretization_->gJEnd(); ++j) {
    // G's boundary condition is derived from the Neumann BC for p:
    // g(0, j) = v(0, j)
    g.at(discretization_->gIBegin() - 1, j) =
        v.at(discretization_->gIBegin() - 1, j);
  }
}

void BoundaryManager::applyRightBoundaryG() {
  FieldVariable &g = discretization_->g();
  const FieldVariable &v = discretization_->v();

  // We start and end in the indices of the ghost cells to
  // prioritize the left-right boundary conditions over the corners
  for (int j = discretization_->gJBegin(); j <= discretization_->gJEnd(); ++j) {
    // G's boundary condition is derived from the Neumann BC for p:
    // g(imax + 1, j) = v(imax + 1, j)
    g.at(discretization_->gIEnd() + 1, j) =
        v.at(discretization_->gIEnd() + 1, j);
  }
}
