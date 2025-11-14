#include "field_variable.h"

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
  auto position_on_grid_x = physicalXPos / _cellSize[0] - _offset[0];
  auto position_on_grid_y = physicalYPos / _cellSize[1] - _offset[1];

  // get the cell index by simply casting to int
  size_t cell_i = position_on_grid_x;
  size_t cell_j = position_on_grid_y;

  // ! panic if the cell is outside the grid
  assert(cell_i < 0);
  assert(cell_i >= _size[0]);
  assert(cell_j < 0);
  assert(cell_j >= _size[1]);

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
