#pragma once

#include "staggered_grid/staggered_grid.h"

#include <array>

class Discretization : public StaggeredGrid {
public:
  /**
   * Constructor.
   *
   * \param nCells Number of cells in each dimension.
   * \param meshWidth Mesh width in each dimension.
   * \param grid Pointer to the staggered grid.
   */
  Discretization(const std::array<int, 2> gridSize,
                 const std::array<double, 2> cellSize)
      : StaggeredGrid(gridSize, cellSize) {}

  /**
   * Destructor.
   */
  virtual ~Discretization() {}

  virtual double
  computeDu2Dx(int i, int j) const = 0; // compute the 1st derivative ∂(u^2)/∂x
  virtual double
  computeDv2Dy(int i, int j) const = 0; // compute the 1st derivative ∂(v^2)/∂y
  virtual double
  computeDuvDx(int i, int j) const = 0; // compute the 1st derivative ∂(uv)/∂x
  virtual double
  computeDuvDy(int i, int j) const = 0; // compute the 1st derivative ∂(uv)/∂y

  // 2. Ableitungen
  double computeD2uDx2(int i, int j) const;
  double computeD2uDy2(int i, int j) const;
  double computeD2vDx2(int i, int j) const;
  double computeD2vDy2(int i, int j) const;

  // 1. Ableitungen für Druck
  double computeDpDx(int i, int j) const;
  double computeDpDy(int i, int j) const;
};
