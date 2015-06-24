#ifndef QUADRATURE2D_HPP
#define QUADRATURE2D_HPP
#include <vector>
#include "vecmath.h"
class Quadrature2D
{
public:
  std::vector<Vector2f> x;
  std::vector<float> w;
  Quadrature2D();
  
  static const Quadrature2D Gauss2;
  static const Quadrature2D Uniform4;
};

#endif