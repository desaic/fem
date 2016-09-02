#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP
#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"
class Quadrature
{
public:
  std::vector<Vector3S> x;
  std::vector<cfgScalar> w;
  Quadrature();
  
  static const Quadrature Gauss2_2D;
  static const Quadrature Gauss4_2D;

  static const Quadrature Gauss2;
  static const Quadrature Uniform4;
};

#endif