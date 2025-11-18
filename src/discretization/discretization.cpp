#include "discretization.h"

 // compute 2nd derivative ∂^2(u)/∂x^2 
double Discretization::computeD2uDx2(int i, int j) const
{
    const double dx = cellSize()[0];
    return (1./(dx*dx)) * (u().at(i+1,j) - 2.*u().at(i,j) + u().at(i-1,j));
}

 // compute 2nd derivative ∂^2(u)/∂y^2 
double Discretization::computeD2uDy2(int i, int j) const
{
    const double dy = cellSize()[1];
    return (1./(dy*dy)) * (u().at(i,j+1) - 2.*u().at(i,j) + u().at(i,j-1));
}

 // compute 2nd derivative ∂^2(v)/∂x^2
double Discretization::computeD2vDx2(int i, int j) const
{
    const double dx = cellSize()[0];
    return (1./(dx*dx)) * (v().at(i+1,j) - 2.*v().at(i,j) + v().at(i-1,j));
}

 // compute 2nd derivative ∂^2(v)/∂y^2
double Discretization::computeD2vDy2(int i, int j) const
{
    const double dy = cellSize()[1];
    return (1./(dy*dy)) * (v().at(i,j+1) - 2.0*v().at(i,j) + v().at(i,j-1));
}

 // compute 1st derivative ∂p/∂x
double Discretization::computeDpDx(int i, int j) const
{
    const double dx = cellSize()[0];
    return (1./dx) * (p().at(i+1,j) - p().at(i,j));
}

 // compute 1st derivative ∂p/∂y
double Discretization::computeDpDy(int i, int j) const
{
    const double dy = cellSize()[1];
    return (1./dy) * (p().at(i,j+1) - p().at(i,j));
}

