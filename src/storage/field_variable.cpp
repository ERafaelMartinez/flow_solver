#include "field_variable.h"

#include <algorithm>
#include <cassert>

// from https://stackoverflow.com/a/4353537
// TODO: Put into its own utility
float lerp(float a, float b, float f) { return a * (1.0 - f) + (b * f); }

FieldVariable::FieldVariable(std::array<int, 2> gridSize,
                             std::array<double, 2> offset,
                             std::array<double, 2> cellSize)
    : Array2D(gridSize), _cellSize(cellSize), _offset(offset) {}

double FieldVariable::interpolateAt(double physicalXPos,
                                    double physicalYPos) const {
  // x and y are being given in phyisical space, so we need to first
  // convert to the corresponding position on our grid
  // TODO: How do we use the positionInCell here?
  auto position_on_grid_x = (physicalXPos / _cellSize[0]) + _offset[0];
  auto position_on_grid_y = (physicalYPos / _cellSize[1]) + _offset[1];

  // get the cell index by simply casting to int
  int cell_i = static_cast<int>(position_on_grid_x);
  int cell_j = static_cast<int>(position_on_grid_y);

  // clamp to valid range
  if (cell_i >= _size[0] - 1) {
    cell_i--;
  }

  if (cell_j >= _size[1] - 1) {
    cell_j--;
  }

  // get the relative position 'inside' the cell
  auto dx = position_on_grid_x - cell_i;
  auto dy = position_on_grid_y - cell_j;

  // interpolate the values from edges of the cell
  // to get the 'exact' value at the 'continuous' point
  auto lowerEdge = lerp(at(cell_i, cell_j), at(cell_i + 1, cell_j), dx);
  auto upperEdge = lerp(at(cell_i, cell_j + 1), at(cell_i + 1, cell_j + 1), dx);
  auto vertical = lerp(lowerEdge, upperEdge, dy);

  return vertical;
}

void FieldVariable::setToZero() {
  for (int i = 0; i < size()[0] * size()[1]; ++i) {
    data()[i] = 0.0;
  }
}

double FieldVariable::maxMagnitude() {
  double max_val = 0.0;
  for (int i = 0; i < size()[0] * size()[1]; ++i) {
    double ith_val = data()[i];
    if (std::fabs(ith_val) > max_val) {
      max_val = std::fabs(ith_val);
    }
  }
  return max_val;
}

// Boundary getters - extract boundary data for MPI communication
std::vector<double> FieldVariable::getLeftBoundary() const {
  std::vector<double> leftCells;
  leftCells.reserve(_size[1]);
  for (int j = 0; j < _size[1]; ++j) {
    leftCells.push_back(at(0, j));
  }
  return leftCells;
}

std::vector<double> FieldVariable::getRightBoundary() const {
  std::vector<double> rightCells;
  rightCells.reserve(_size[1]);
  for (int j = 0; j < _size[1]; ++j) {
    rightCells.push_back(at(_size[0] - 1, j));
  }
  return rightCells;
}

std::vector<double> FieldVariable::getTopBoundary() const {
  std::vector<double> topCells;
  topCells.reserve(_size[0]);
  for (int i = 0; i < _size[0]; ++i) {
    topCells.push_back(at(i, _size[1] - 1));
  }
  return topCells;
}

std::vector<double> FieldVariable::getBottomBoundary() const {
  std::vector<double> bottomCells;
  bottomCells.reserve(_size[0]);
  for (int i = 0; i < _size[0]; ++i) {
    bottomCells.push_back(at(i, 0));
  }
  return bottomCells;
}

// Boundary setters - update boundary data after MPI communication
void FieldVariable::setLeftBoundary(const std::vector<double> &data) {
  assert(data.size() == _size[1]);
  for (int j = 0; j < _size[1]; ++j) {
    at(0, j) = data[j];
  }
}

void FieldVariable::setRightBoundary(const std::vector<double> &data) {
  assert(data.size() == _size[1]);
  for (int j = 0; j < _size[1]; ++j) {
    at(_size[0] - 1, j) = data[j];
  }
}

void FieldVariable::setTopBoundary(const std::vector<double> &data) {
  assert(data.size() == _size[0]);
  for (int i = 0; i < _size[0]; ++i) {
    at(i, _size[1] - 1) = data[i];
  }
}

void FieldVariable::setBottomBoundary(const std::vector<double> &data) {
  assert(data.size() == _size[0]);
  for (int i = 0; i < _size[0]; ++i) {
    at(i, 0) = data[i];
  }
}
