#include "array2d.h"

#include <cassert>

Array2D::Array2D(std::array<int, 2> size) : _size(size)
{
    // allocate data, initialize to 0
    _data.resize(_size[0] * _size[1], 0.0);
}

//! get the size
std::array<int, 2> Array2D::size() const
{
    return _size;
}

double &Array2D::at(int i, int j)
{
    const int index = j * _size[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < _size[0]);
    assert(0 <= j && j < _size[1]);
    assert(j * _size[0] + i < (int)_data.size());

    return _data[index];
}

double Array2D::at(int i, int j) const
{
    const int index = j * _size[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < _size[0]);
    assert(0 <= j && j < _size[1]);
    assert(j * _size[0] + i < (int)_data.size());

    return _data[index];
}
