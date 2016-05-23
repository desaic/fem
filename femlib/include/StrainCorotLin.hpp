#ifndef STRAINCOROTLIN_HPP
#define STRAINCOROTLIN_HPP
#include "StrainLin.hpp"

#include <Eigen/Dense>

class StrainCorotLin :public StrainLin{
public:
  StrainCorotLin();
  virtual cfgScalar getEnergy(const Matrix3S & F);
  virtual Matrix3S getPK1(const Matrix3S & F);
  virtual Matrix3S getdPdx(const Matrix3S & F, const Matrix3S & dF);
  MatrixXS EMatrix();
};

#endif
