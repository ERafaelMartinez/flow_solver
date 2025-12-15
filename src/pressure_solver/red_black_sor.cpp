#include "red_black_sor.h"

/* This file implements the RedBlackSORPressureSolver class, which is a specific
   implementation of the PressureSolver interface. It uses a modified version of
   the SOR iterative method to solve the pressure Poisson equation in a stable
   manner that works in parallel scenarios
*/

// Constructor
RedBlackSORPressureSolver::RedBlackSORPressureSolver(
    std::shared_ptr<Discretization> discretization,
    std::shared_ptr<Partitioning> partitioning, double convergence_tol,
    int max_iterations, double omega)
    : PressureSolver(discretization, partitioning, convergence_tol,
                     max_iterations),
      omega_(omega) {}

// Calculate pressure for one iteration
void RedBlackSORPressureSolver::calcPressureIter() {
  auto cellSize = discretization_->cellSize();
  const double dx = cellSize[0];
  const double dy = cellSize[1];
  const double idx2 = 1.0 / (dx * dx);
  const double idy2 = 1.0 / (dy * dy);

  // Precompute coefficient
  const double coeff = dx * dx * dy * dy / (2.0 * (dx * dx + dy * dy));

  int redStartingIndex = 0;
  int blackStartingIndex = 1;

  // 1. First do the 'red' cells
  doColorIteration(GridColor::RED, coeff, idx2, idy2);

  // 2. Then exchange the values
  // It should be fine to exchange all values of p here.
  dataExchanger_->exchange(discretization_->p());

  // 3. Then do the 'black' cells
  doColorIteration(GridColor::BLACK, coeff, idx2, idy2);

  // 4. Then exchange the values again
  // No need to explicitely exchange values here since this will be
  // automatically done in the abstract.
}

void RedBlackSORPressureSolver::doColorIteration(
    RedBlackSORPressureSolver::GridColor color, const double coeff,
    const double idx2, const double idy2) {

  // Since the grid follows a checkerboard pattern, we need to shift the columns
  // on each row to follow this pattern.
  //    0 1 2 3 4
  //    ----------
  // 3| ⚫️🔴⚫️🔴⚫️
  // 2| 🔴⚫️🔴⚫️🔴
  // 1| ⚫️🔴⚫️🔴⚫️
  // 0| 🔴⚫️🔴⚫️🔴
  int colShiftOnEvenRows = 0, colShiftOnOddRows = 0;
  if (color == GridColor::RED) {
    colShiftOnEvenRows = 0;
    colShiftOnOddRows = 1;
  } else if (color == GridColor::BLACK) {
    colShiftOnEvenRows = 1;
    colShiftOnOddRows = 0;
  } else {
    throw std::invalid_argument("Unknown color");
  }

  FieldVariable &p = discretization_->p();
  const FieldVariable &rhs = discretization_->rhs();

  for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
    int colOffset = j % 2 == 0 ? colShiftOnEvenRows : colShiftOnOddRows;

    for (int i = discretization_->pIBegin() + colOffset;
         i <= discretization_->pIEnd(); i += 2) {
      // Compute the Gauss-Seidel target value
      const double p_GS =
          coeff * ((p.at(i - 1, j) + p.at(i + 1, j)) * idx2 +
                   (p.at(i, j - 1) + p.at(i, j + 1)) * idy2 - rhs.at(i, j));

      const double p_old = p.at(i, j);
      const double p_new = (1.0 - omega_) * p_old + omega_ * p_GS;
      p.at(i, j) = p_new;
    }
  }
}
