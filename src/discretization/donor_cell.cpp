#include "donor_cell.h"
#include <cassert>

DonorCell::DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha)
    : Discretization(nCells, meshWidth), alpha_(alpha)
{
    assert(alpha_ >= 0.0 && alpha_ <= 1.0 && "alpha must be in [0,1]");
}

// compute the 1st derivative ∂(u^2)/∂x
double DonorCell::computeDu2Dx(int i, int j) const
{
    const double dx = meshWidth_[0];
    const double uRight = (u(i,j) + u(i+1,j)) / 2.;
    const double uLeft = (u(i-1,j) + u(i,j)) / 2.;
    const double uDiffRight = (u(i,j) - u(i+1,j)) / 2.;
    const double uDiffLeft = (u(i-1,j) - u(i,j)) / 2.;
    CentralDifferences cent(nCells_, meshWidth_);
    const double central = cent.computeDu2Dx(i,j);
    return central + alpha_ * 1./dx * ((abs(uRight) * uDiffRight) - (abs(uLeft) * uDiffLeft) );
}
// compute the 1st derivative ∂(v^2)/∂y
double DonorCell::computeDv2Dy(int i, int j) const 
{
    const double dy = meshWidth_[1];
    const double vTop = (v(i,j) + v(i,j+1)) / 2.;
    const double vBottom = (v(i,j-1) + v(i,j)) / 2.;
    const double vDiffTop = (v(i,j) - v(i,j+1)) / 2.;
    const double vDiffBottom = (v(i,j-1) - v(i,j)) / 2.;
    CentralDifferences cent(nCells_, meshWidth_);
    const double central = cent.computeDv2Dy(i,j);
    
    return central + alpha_ * 1./dy * ((abs(vTop) * vDiffTop) - (abs(vBottom) * vDiffBottom) );
}
// compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const
{
    const double dx = meshWidth_[0];
    const double uTopRight = (u(i,j) + u(i,j+1)) / 2.;
    const double uTopLeft = (u(i-1,j) + u(i-1,j+1)) / 2.;
    const double vDiffRight = (v(i,j) - v(i+1,j)) / 2.;
    const double vDiffLeft = (v(i-1,j) - v(i,j)) / 2.;
    CentralDifferences cent(nCells_, meshWidth_);
    const double central = cent.computeDuvDx(i,j);
    return central + alpha_ * 1./dx * ((abs(uTopRight) * vDiffRight) - (abs(uTopLeft) * vDiffLeft));
}
// compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const 
{
    const double dy = meshWidth_[1];
    const double vTopRight = (v(i,j) + v(i+1,j)) / 2.;
    const double vBottomRight = (v(i,j-1) + v(i+1,j-1)) / 2.;
    const double uDiffTop = (u(i,j) - u(i,j+1)) / 2.;
    const double uDiffBottom = (u(i,j-1) - u(i,j)) / 2.;
    CentralDifferences cent(nCells_, meshWidth_);
    const double central = cent.computeDuvDy(i,j);
    return central + alpha_ * 1./dy * ((abs(vTopRight) * uDiffTop) - (abs(vBottomRight) * uDiffBottom) );
}
