#include "donor_cell.h"
#include "central_differences.h"

#include <cassert>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth,
                     double alpha)
    : Discretization(nCells, meshWidth), alpha_(alpha) {
  assert(alpha_ >= 0.0 && alpha_ <= 1.0 && "alpha must be in [0,1]");
}

// compute the 1st derivative ∂(u^2)/∂x
double DonorCell::computeDu2Dx(int i, int j) const {
  const double dx = cellSize()[0];
  const double uRight = (u().at(i, j) + u().at(i + 1, j)) / 2.;
  const double uLeft = (u().at(i - 1, j) + u().at(i, j)) / 2.;
  const double uDiffRight = (u().at(i, j) - u().at(i + 1, j)) / 2.;
  const double uDiffLeft = (u().at(i - 1, j) - u().at(i, j)) / 2.;

  // TODO: remove weird dep
  CentralDifferences cent(gridSize(), cellSize());
  const double central = cent.computeDu2Dx(i, j);
  const double correction =
      ((abs(uRight) * uDiffRight) - (abs(uLeft) * uDiffLeft));
  return central + (1. / dx) * alpha_ * correction;
}
// compute the 1st derivative ∂(v^2)/∂y
double DonorCell::computeDv2Dy(int i, int j) const {
  const double dy = cellSize()[1];
  const double vTop = (v().at(i, j) + v().at(i, j + 1)) / 2.;
  const double vBottom = (v().at(i, j - 1) + v().at(i, j)) / 2.;
  const double vDiffTop = (v().at(i, j) - v().at(i, j + 1)) / 2.;
  const double vDiffBottom = (v().at(i, j - 1) - v().at(i, j)) / 2.;

  // TODO: remove weird dep
  CentralDifferences cent(gridSize(), cellSize());
  const double central = cent.computeDv2Dy(i, j);
  const double correction =
      ((abs(vTop) * vDiffTop) - (abs(vBottom) * vDiffBottom));

  return central + (1. / dy) * alpha_ * correction;
}
// compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const {
  const double dx = cellSize()[0];
  const double uTopRight = (u().at(i, j) + u().at(i, j + 1)) / 2.;
  const double uTopLeft = (u().at(i - 1, j) + u().at(i - 1, j + 1)) / 2.;
  const double vDiffRight = (v().at(i, j) - v().at(i + 1, j)) / 2.;
  const double vDiffLeft = (v().at(i - 1, j) - v().at(i, j)) / 2.;

  // TODO: remove weird dep
  CentralDifferences cent(gridSize(), cellSize());
  const double central = cent.computeDuvDx(i, j);
  const double correction =
      (abs(uTopRight) * vDiffRight) - (abs(uTopLeft) * vDiffLeft);

  return central + (1. / dx) * alpha_ * correction;
}
// compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const {
  const double dy = cellSize()[1];
  const double vTopRight = (v().at(i, j) + v().at(i + 1, j)) / 2.;
  const double vBottomRight = (v().at(i, j - 1) + v().at(i + 1, j - 1)) / 2.;
  const double uDiffTop = (u().at(i, j) - u().at(i, j + 1)) / 2.;
  const double uDiffBottom = (u().at(i, j - 1) - u().at(i, j)) / 2.;

  // TODO: remove weird dep
  CentralDifferences cent(gridSize(), cellSize());
  const double central = cent.computeDuvDy(i, j);
  const double correction =
      (abs(vTopRight) * uDiffTop) - (abs(vBottomRight) * uDiffBottom);

  return central + (1. / dy) * alpha_ * correction;
}
