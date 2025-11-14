#pragma once

#include <array>
#include "discretization.h"  

class CentralDifferences : public Discretization 
{
public:
/**
 * Constructor.
 * 
 * \param nCells Number of cells in each dimension.
 * \param meshWidth Mesh width in each dimension.
 */
CentralDifferences(std::array<int,2> nCells, std::array<double,2> meshWidth);

double computeDu2Dx(int i, int j) const override;
double computeDv2Dy(int i, int j) const override;
double computeDuvDx(int i, int j) const override;
double computeDuvDy(int i, int j) const override;
};
