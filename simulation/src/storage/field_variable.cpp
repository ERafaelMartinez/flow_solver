#include "field_variable.h"

#include <algorithm>
#include <cassert>
#include <cmath>

//! Linear interpolation between two values
// x1, x2: values to interpolate between
// x: relative position between x1 and x2 (0 <= x <= 1)
double linear_interpolation(double x1, double x2, double x) {
  return x1 * (1.0 - x) + (x2 * x);
}

FieldVariable::FieldVariable(std::array<int, 2> gridSize,
                             std::array<double, 2> origin,
                             std::array<double, 2> cellSize)
    : Array2D(gridSize), _cellSize(cellSize), _origin(origin) {}

double FieldVariable::interpolateAt(double physicalXPos,
                                    double physicalYPos) const {
  // x and y are given in the reference frame of the physical domain
  // they must be converted to the reference frame of the grid (shifted by
  // _origin)

  // first scale the physical coordinates to relative grid coordinates
  double position_on_domain_grid_x =
      (physicalXPos / _cellSize[0]); // units: {cells}
  double position_on_domain_grid_y =
      (physicalYPos / _cellSize[1]); // units: {cells}
  // apply the origin offset, which is given in cell units
  position_on_domain_grid_x -= _origin[0];
  position_on_domain_grid_y -= _origin[1];

  // obtain the cell index by taking the floor of the value
  int cell_i = static_cast<int>(position_on_domain_grid_x);
  int cell_j = static_cast<int>(position_on_domain_grid_y);

  // if the point is out of the simulation domain, this will cause
  // 'invalid' cell indices. The valid domain indices range over
  // [0, _array_size - 2] --> this is because of the ghost cells
  // used.
  // In case the point is out of the domain,
  // we assign the values of the nearest valid cell
  if (cell_i > _size[0] - 2) {
    cell_i = _size[0] - 2;
  }
  if (cell_i < 0) {
    cell_i = 0;
  }
  if (cell_j > _size[1] - 2) {
    cell_j = _size[1] - 2;
  }
  if (cell_j < 0) {
    cell_j = 0;
  }

  // get the relative position 'inside' the cell
  double dx = position_on_domain_grid_x - cell_i;
  double dy = position_on_domain_grid_y - cell_j;

  // interpolate the values from edges of the cell linearly

  // horizontal interpolation:
  double varAtXLowerEdge =
      linear_interpolation(at(cell_i, cell_j), at(cell_i + 1, cell_j), dx);
  double varAtXUpperEdge = linear_interpolation(at(cell_i, cell_j + 1),
                                                at(cell_i + 1, cell_j + 1), dx);
  double varAtXYPoint =
      linear_interpolation(varAtXLowerEdge, varAtXUpperEdge, dy);
  return varAtXYPoint;
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
std::vector<double>
FieldVariable::getColumnValues(int columnIndex,
                               std::array<int, 2> rowRange) const {
  std::vector<double> columnValues;
  columnValues.reserve(_size[1]);
  for (int j = rowRange[0]; j < rowRange[1]; ++j) {
    columnValues.push_back(at(columnIndex, j));
  }
  return columnValues;
}

std::vector<double>
FieldVariable::getRowValues(int rowIndex,
                            std::array<int, 2> columnRange) const {
  std::vector<double> rowValues;
  rowValues.reserve(_size[0]);
  for (int i = columnRange[0]; i < columnRange[1]; ++i) {
    rowValues.push_back(at(i, rowIndex));
  }
  return rowValues;
}

// Boundary setters - update boundary data after MPI communication
void FieldVariable::setColumnValues(int columnIndex,
                                    std::array<int, 2> rowRange,
                                    const std::vector<double> &columnValues) {
  // assert that the number of values matches the number of rows
  assert(columnValues.size() == rowRange[1] - rowRange[0]);
  // assert that the row range is valid
  assert(rowRange[0] >= 0 && rowRange[1] <= _size[1]);

  for (int j = rowRange[0]; j < rowRange[1]; ++j) {
    at(columnIndex, j) = columnValues[j - rowRange[0]];
  }
}

void FieldVariable::setRowValues(int rowIndex, std::array<int, 2> columnRange,
                                 const std::vector<double> &rowValues) {
  // assert that the number of values matches the number of columns
  assert(rowValues.size() == columnRange[1] - columnRange[0]);
  // assert that the column range is valid
  assert(columnRange[0] >= 0 && columnRange[1] <= _size[0]);

  for (int i = columnRange[0]; i < columnRange[1]; ++i) {
    at(i, rowIndex) = rowValues[i - columnRange[0]];
  }
}
