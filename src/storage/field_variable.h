#pragma once

#include "array2d.h"

#include <array>

class FieldVariable : public Array2D
{
public:
    FieldVariable(std::array<int, 2> gridSize,
                  std::array<double, 2> offset,
                  std::array<double, 2> cellSize);

    double interpolateAt(double x, double y) const;

    double maxMagnitude();

private:
    const std::array<double, 2> _cellSize;
    const std::array<double, 2> _offset;
};
