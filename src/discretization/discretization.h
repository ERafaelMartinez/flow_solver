#pragma once

#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include <array> 
#include "staggered_grid.h"

class Discretization : public StaggeredGrid 
{
protected:
    std::array<int, 2> nCells_;     // Number of cells in each dimension
    std::array<double, 2> meshWidth_; // dx, dy

public:
/**
 * Constructor.
 * 
 * \param nCells Number of cells in each dimension.
 * \param meshWidth Mesh width in each dimension.
 * \param grid Pointer to the staggered grid.
 */
Discretization(const std::array<int,2>& nCells, 
               const std::array<double,2>& meshWidth);

/**
 * Destructor.
 */
virtual ~Discretization() {}


virtual double computeDu2Dx(int i, int j) const = 0;   // compute the 1st derivative ∂(u^2)/∂x
virtual double computeDv2Dy(int i, int j) const = 0;   // compute the 1st derivative ∂(v^2)/∂y 
virtual double computeDuvDx(int i, int j) const = 0;   // compute the 1st derivative ∂(uv)/∂x
virtual double computeDuvDy(int i, int j) const = 0;   // compute the 1st derivative ∂(uv)/∂y 

// 2. Ableitungen
double computeD2uDx2(int i, int j) const;
double computeD2uDy2(int i, int j) const;
double computeD2vDx2(int i, int j) const;
double computeD2vDy2(int i, int j) const;

// 1. Ableitungen für Druck
double computeDpDx(int i, int j) const;
double computeDpDy(int i, int j) const;

};

#endif
