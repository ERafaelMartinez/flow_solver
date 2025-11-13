#include "staggered_grid.h"

StaggeredGrid::StaggeredGrid(std::array<int, 2> gridSize,
                             std::array<double, 2> cellSize)
    : _gridSize(gridSize), _cellSize(cellSize),
      _u(getVarGridSize(gridSize, _uBoundaries), {0.0, 0.5}, cellSize),
      _v(getVarGridSize(gridSize, _vBoundaries), {0.5, 1.0}, cellSize),
      _p(getVarGridSize(gridSize, _pBoundaries), {0.5, 0.5}, cellSize),
      _rhs(getVarGridSize(gridSize, _rhsBoundaries), {0.5, 0.5}, cellSize),
      _g(getVarGridSize(gridSize, _gBoundaries), {0.0, 0.0}, cellSize),
      _f(getVarGridSize(gridSize, _fBoundaries), {0.5, 1.0}, cellSize) {}
