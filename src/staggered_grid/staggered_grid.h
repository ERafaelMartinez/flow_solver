#include "field_variable.h"

enum INDEX_TYPE { I_BEGIN, I_END, J_BEGIN, J_END };

class StaggeredGrid {
private:
  // {top, right, bottom, left}
  std::array<int, 4> _uBoundaries = {1, 0, 1, 1};
  std::array<int, 4> _vBoundaries = {0, 1, 1, 1};
  std::array<int, 4> _pBoundaries = {1, 1, 1, 1};
  std::array<int, 4> _rhsBoundaries = {1, 1, 1, 1};
  std::array<int, 4> _fBoundaries = {1, 0, 1, 1};
  std::array<int, 4> _gBoundaries = {0, 1, 1, 1};

  FieldVariable _u;
  FieldVariable _v;
  FieldVariable _p;
  FieldVariable _rhs;
  FieldVariable _g;
  FieldVariable _f;

  std::array<int, 2> _gridSize;
  std::array<double, 2> _cellSize;

  std::array<int, 2> getVarGridSize(std::array<int, 2> gridSize,
                                    std::array<int, 4> boundaries) const {
    return {gridSize[0] + _uBoundaries[1] + _uBoundaries[3],
            gridSize[1] + _uBoundaries[0] + _uBoundaries[2]};
  }

  int getInnerIndex(FieldVariable var, std::array<int, 4> boundaries,
                    INDEX_TYPE index_type) const {
    switch (index_type) {
    case I_BEGIN: {
      return var.size()[0] + boundaries[3];
    }
    case I_END: {
      return var.size()[0] + boundaries[1];
    }
    case J_BEGIN: {
      return var.size()[1] + boundaries[2];
    }
    case J_END: {
      return var.size()[1] + boundaries[0];
    }
    }
  }

public:
  StaggeredGrid(std::array<int, 2> gridSize, std::array<double, 2> cellSize);

  std::array<int, 2> gridSize() const { return _gridSize; }
  std::array<double, 2> cellSize() const { return _cellSize; }

  FieldVariable u() const { return _u; }
  FieldVariable v() const { return _v; }
  FieldVariable p() const { return _p; }
  FieldVariable rhs() const { return _rhs; }
  FieldVariable g() const { return _g; }
  FieldVariable f() const { return _f; }

  int uIBegin() const { return getInnerIndex(_u, _uBoundaries, I_BEGIN); }
  int uIEnd() const { return getInnerIndex(_u, _uBoundaries, I_END); }
  int uJBegin() const { return getInnerIndex(_u, _uBoundaries, J_BEGIN); }
  int uJEnd() const { return getInnerIndex(_u, _uBoundaries, J_END); }

  int vIBegin() const { return getInnerIndex(_v, _vBoundaries, I_BEGIN); }
  int vIEnd() const { return getInnerIndex(_v, _vBoundaries, I_END); }
  int vJBegin() const { return getInnerIndex(_v, _vBoundaries, J_BEGIN); }
  int vJEnd() const { return getInnerIndex(_v, _vBoundaries, J_END); }

  int pIBegin() const { return getInnerIndex(_p, _pBoundaries, I_BEGIN); }
  int pIEnd() const { return getInnerIndex(_p, _pBoundaries, I_END); }
  int pJBegin() const { return getInnerIndex(_p, _pBoundaries, J_BEGIN); }
  int pJEnd() const { return getInnerIndex(_p, _pBoundaries, J_END); }

  int rhsIBegin() const { return getInnerIndex(_rhs, _rhsBoundaries, I_BEGIN); }
  int rhsIEnd() const { return getInnerIndex(_rhs, _rhsBoundaries, I_END); }
  int rhsJBegin() const { return getInnerIndex(_rhs, _rhsBoundaries, J_BEGIN); }
  int rhsJEnd() const { return getInnerIndex(_rhs, _rhsBoundaries, J_END); }

  int gIBegin() const { return getInnerIndex(_g, _gBoundaries, I_BEGIN); }
  int gIEnd() const { return getInnerIndex(_g, _gBoundaries, I_END); }
  int gJBegin() const { return getInnerIndex(_g, _gBoundaries, J_BEGIN); }
  int gJEnd() const { return getInnerIndex(_g, _gBoundaries, J_END); }

  int fIBegin() const { return getInnerIndex(_f, _fBoundaries, I_BEGIN); }
  int fIEnd() const { return getInnerIndex(_f, _fBoundaries, I_END); }
  int fJBegin() const { return getInnerIndex(_f, _fBoundaries, J_BEGIN); }
  int fJEnd() const { return getInnerIndex(_f, _fBoundaries, J_END); }
};
