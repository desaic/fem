#ifndef QUADRATURE2D_HPP
#define QUADRATURE2D_HPP
#include <vector>
#include "cfgDefs.h"

class Quadrature2D
{
public:
  std::vector<Vector2S> x;
  std::vector<cfgScalar> w;
  Quadrature2D();
  
  static const Quadrature2D Gauss2;
  static const Quadrature2D Uniform4;
};

#endif