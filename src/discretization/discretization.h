#include "staggered_grid/staggered_grid.h"

class Discretization : public StaggeredGrid {
  Discretization(std::array<int, 2> gridSize, std::array<double, 2> cellSize)
      : StaggeredGrid(gridSize, cellSize) {}
};
