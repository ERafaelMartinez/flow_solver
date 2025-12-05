#include "field_variable.h"

#include <algorithm>
#include <cassert>
#include <cmath>

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

// row/column getters - extract boundary data for MPI communication
std::vector<double> FieldVariable::getColumnValues(int columnIndex, std::array<int, 2> rowRange) const {
  std::vector<double> columnValues;
  columnValues.reserve(_size[1]);
  for (int j = rowRange[0]; j < rowRange[1]; ++j) {
    columnValues.push_back(at(columnIndex, j));
  }
  return columnValues;
}

std::vector<double> FieldVariable::getRowValues(int rowIndex, std::array<int, 2> columnRange) const {
  std::vector<double> rowValues;
  rowValues.reserve(_size[0]);
  for (int i = columnRange[0]; i < columnRange[1]; ++i) {
    rowValues.push_back(at(i, rowIndex));
  }
  return rowValues;
}

// Boundary setters - update boundary data after MPI communication
void FieldVariable::setColumnValues(int columnIndex, std::array<int, 2> rowRange, const std::vector<double> &columnValues) {
  // assert that the number of values matches the number of rows
  assert(columnValues.size() == rowRange[1] - rowRange[0]);
  // assert that the row range is valid
  assert(rowRange[0] >= 0 && rowRange[1] <= _size[1]);

  for (int j = rowRange[0]; j < rowRange[1]; ++j) {
    at(columnIndex, j) = columnValues[j];
  }
}

void FieldVariable::setRowValues(int rowIndex, std::array<int, 2> columnRange, const std::vector<double> &rowValues) {
  // assert that the number of values matches the number of columns
  assert(rowValues.size() == columnRange[1] - columnRange[0]);
  // assert that the column range is valid
  assert(columnRange[0] >= 0 && columnRange[1] <= _size[0]);

  for (int i = columnRange[0]; i < columnRange[1]; ++i) {
    at(i, rowIndex) = rowValues[i];
  }
}

