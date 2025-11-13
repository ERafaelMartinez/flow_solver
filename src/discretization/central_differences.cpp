#include "central_differences.h"


CentralDifferences::CentralDifferences(std::array<int,2> nCells, std::array<double,2> meshWidth)
    : Discretization(nCells, meshWidth)
{  
}


// compute the 1st derivative ∂(u^2)/∂x 
double CentralDifferences::computeDu2Dx(int i, int j) const
{
    const double dx = meshWidth_[0];
    const double uRight = (u(i+1,j) + u(i,j)) / 2.;
    const double uLeft = (u(i,j) + u(i-1,j)) / 2.;
    return 1./dx * (uRight*uRight - uLeft*uLeft);
}

// compute the 1st derivative ∂(v^2)/∂y 
double CentralDifferences::computeDv2Dy(int i, int j) const 
{
    const double dy = meshWidth_[1];
    const double vTop = (v(i,j+1) + v(i,j)) / 2.;
    const double vBottom = (v(i,j) + v(i,j-1)) / 2.;
    return 1./dy * (vTop*vTop - vBottom*vBottom);
}

// compute the 1st derivative ∂(uv)/∂x 
double CentralDifferences::computeDuvDx(int i, int j) const 
{
    const double dx = meshWidth_[0];     // mesh width in x direction, δx
    const double uTopRight  = (u(i,j) + u(i,j+1)) / 2.;         // u at top right corner of cell
    const double uTopLeft  = (u(i-1,j) + u(i-1,j+1)) / 2.;      // u at top left corner of cell
    const double vTopRight  = (v(i,j) + v(i+1,j)) / 2.;         // v at top right corner of cell
    const double vTopLeft  = (v(i-1,j) + v(i,j)) / 2.;          // v at top left corner of cell
    
    return 1./dx * (uTopRight * vTopRight - uTopLeft * vTopLeft);
}

// compute the 1st derivative ∂(uv)/∂y
double CentralDifferences::computeDuvDy(int i, int j) const
{
    const double dy = meshWidth_[1];      // mesh width in y direction, δy
    
    const double vTopRight  = (v(i,j) + v(i+1,j)) / 2.;         // v at top right corner of cell
    const double vBottomRight = (v(i,j-1) + v(i+1,j-1)) / 2.;   // v at bottom right corner of cell
    const double uTopRight  = (u(i,j) + u(i,j+1)) / 2.;         // u at top right corner of cell
    const double uBottomRight = (u(i,j-1) + u(i,j)) / 2.;       // u at bottom right corner of cell
    
    return 1./dy * (vTopRight * uTopRight - vBottomRight * uBottomRight);
}