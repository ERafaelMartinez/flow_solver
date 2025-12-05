#pragma once

#include "array2d.h"

#include <array>

class FieldVariable : public Array2D {
public:
  FieldVariable(std::array<int, 2> gridSize, std::array<double, 2> offset,
                std::array<double, 2> cellSize);

  double interpolateAt(double x, double y) const;

  double maxMagnitude();

  void setToZero();

  // Boundary accessors for MPI communication
  std::vector<double> getLeftBoundary() const;
  std::vector<double> getRightBoundary() const;
  std::vector<double> getTopBoundary() const;
  std::vector<double> getBottomBoundary() const;

  // Boundary setters for MPI communication
  void setLeftBoundary(const std::vector<double> &data);
  void setRightBoundary(const std::vector<double> &data);
  void setTopBoundary(const std::vector<double> &data);
  void setBottomBoundary(const std::vector<double> &data);

private:
  const std::array<double, 2> _cellSize;
  const std::array<double, 2> _offset;
};
