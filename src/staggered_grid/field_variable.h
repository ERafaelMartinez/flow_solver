#pragma once

#ifndef FIELDVARIABLE_H
#define FIELDVARIABLE_H

#include <array> 
#include "array2d.h"

class FieldVariable : public Array2D
{
public:
FieldVariable (std::array< int, 2 > size, std::array< double, 2 > origin, std::array< double, 2 > meshWidth)
   
double 	interpolateAt (double x, double y) const;   // get the value at the Cartesian coordinate (x,y)

private:
    std::array< double, 2 > origin_;     // Origin of the field variable
    std::array< double, 2 > meshWidth_;  // Mesh width in each dimension
};

#endif