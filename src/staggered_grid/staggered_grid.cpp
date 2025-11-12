#include "staggered_grid.h"

StaggeredGrid::(std::array<int,2> nCells, std::array<double,2> meshWidth);
       : FieldVariable(std::array<int,2>{nCells[0]+2, nCells[1]+2},   // size including ghost cells
            std::array<double,2>{0.0, 0.0},                 // origin at (0,0)
            meshWidth                                      // mesh width
        ),
{
}