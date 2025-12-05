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

  // row/column getters for MPI communication
  std::vector<double> getColumnValues(int columnIndex, std::array<int, 2> rowRange) const;
  std::vector<double> getRowValues(int rowIndex, std::array<int, 2> columnRange) const;

  // bulk row/column setters for MPI communication
  void setColumnValues(int columnIndex, std::array<int, 2> rowRange, const std::vector<double> &columnValues);
  void setRowValues(int rowIndex, std::array<int, 2> columnRange, const std::vector<double> &rowValues);

private:
  const std::array<double, 2> _cellSize;
  const std::array<double, 2> _offset;
};
