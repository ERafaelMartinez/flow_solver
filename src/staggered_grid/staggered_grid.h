#pragma once

#ifndef STAGGEREDGRID_H
#define STAGGEREDGRID_H

#include <array> 
#include "field_variable.h"

class StaggeredGrid : public FieldVariable
{
public:
    StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth);
       
    const std::array< double, 2 > 	meshWidth () const;   // get the mesh width, i.e. the length of a single cell in x and y direction
    const std::array< int, 2 > 	nCells () const;          // get number of cells in each coordinate direction
    const FieldVariable & 	u () const;                   // get a reference to field variable u
    const FieldVariable & 	v () const;                   // get a reference to field variable v
    const FieldVariable & 	p () const;                   // get a reference to field variable p
  
    double u(int i, int j) const;                         // access value of u in element (i,j)
    double& u(int i, int j);                              // access value of u in element (x,y)
    double v(int i, int j) const;                         // access value of v in element (i,j)
    double& v(int i, int j);                              // access value of v in element (x,y)
    double p(int i, int j) const;                         // access value of p in element (i,j)
    double& p(int i, int j);                              // access value of p in element (x,y)
    double & rhs (int i, int j);                          // access value of rhs in element (i,j)
    double & f (int i, int j);                            // access value of f in element (i,j)
    double & g (int i, int j);                            // access value of g in element (i,j)
    double 	dx () const;                                  // get the mesh width in x-direction, δx 
    double 	dy () const;                                  // get the mesh width in y-direction, δy
    int uIBegin () const;                                 // first valid index for u in x direction
    int uIEnd () const;                                   // one after last valid index for u in x direction
    int uJBegin () const;                                 // first valid index for u in y direction
    int uJEnd () const;                                   // one after last valid index for u in y
    int vIBegin () const;                                 // first valid index for v in x direction
    int vIEnd () const;                                   // one after last valid index for v in x direction
    int vJBegin () const;                                 // first valid index for v in y direction
    int vJEnd () const;                                   // one after last valid index for v in y
    int pIBegin () const;                                 // first valid index for p in x direction
    int pIEnd () const;                                   // one after last valid index for p in x direction
    int pJBegin () const;                                 // first valid index for p in y direction
    int pJEnd () const;                                   // one after last valid index for p in y


protected:
    std::array<int,2> nCells_;
    std::array<double,2> meshWidth_;
    FieldVariable u_, v_, p_, rhs_, f_, g_;
};

#endif  // STAGGEREDGRID_H
