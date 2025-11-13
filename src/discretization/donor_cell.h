#pragma once

#ifndef DONORCELL_H
#define DONORCELL_H

#include <cmath>
#include <array>
#include "discretization.h"  
#include "staggered_grid.h"
#include "central_differences.h"

class DonorCell : public Discretization {
    double alpha_; // blending parameter
public:
/**
 * Constructor.
 * 
 * \param nCells Number of cells in each dimension.
 * \param meshWidth Mesh width in each dimension.
 * \param alpha Blending parameter.
 */
DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha);

double computeDu2Dx(int i, int j) const override;
double computeDv2Dy(int i, int j) const override;
double computeDuvDx(int i, int j) const override;
double computeDuvDy(int i, int j) const override;
};

#endif
