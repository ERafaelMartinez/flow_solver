#include "staggered_grid/array2d.h"

#include <cassert>

Array2D::Array2D(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0]*size_[1], 0.0);
}

//! get the size
std::array<int,2> Array2D::size() const
{
  return size_;
}

double &Array2D::operator()(int i, int j)
{
  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  const int index = j*size_[0] + i;

  return data_[index];
}

double Array2D::operator()(int i, int j) const
{
  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());
  
  const int index = j*size_[0] + i;

  return data_[index];
}