#include "field_variable.h"

#include <algorithm>
#include <cassert>

// from https://stackoverflow.com/a/4353537
// TODO: Put into its own utility
float lerp(float a, float b, float f) { return a * (1.0 - f) + (b * f); }

FieldVariable::FieldVariable(std::array<int, 2> gridSize,
                             std::array<double, 2> offset,
                             std::array<double, 2> cellSize)
    : Array2D(gridSize), _offset(offset), _cellSize(cellSize) {}

double FieldVariable::interpolateAt(double physicalXPos,
                                    double physicalYPos) const {
  // x and y are being given in phyisical space, so we need to first
  // convert to the corresponding position on our grid
  // TODO: How do we use the positionInCell here?
  auto position_on_grid_x = (physicalXPos / _cellSize[0]) + _offset[0];
  auto position_on_grid_y = (physicalYPos / _cellSize[1]) + _offset[1];

  // get the cell index by simply casting to int
  std::size_t cell_i = position_on_grid_x;
  std::size_t cell_j = position_on_grid_y;

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


double FieldVariable::maxMagnitude() {
  double max_val = 0.0;
  for (int i = 0; i <= size()[0]*size()[1]; ++i) {
    double ith_val = data()[i];
    if (std::abs(ith_val) > max_val) {max_val = std::abs(ith_val);}
  }
  return max_val;
}
