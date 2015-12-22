#ifndef STRAINCOROTLIN_HPP
#define STRAINCOROTLIN_HPP
#include "StrainLin.hpp"

#include <Eigen/Dense>

class StrainCorotLin :public StrainLin{
public:
  StrainCorotLin();
  virtual float getEnergy(const Eigen::Matrix3f & F);
  virtual Eigen::Matrix3f getPK1(const Eigen::Matrix3f & F);
  virtual Eigen::Matrix3f getdPdx(const Eigen::Matrix3f & F, const Eigen::Matrix3f & dF);
  Eigen::MatrixXf EMatrix();
};

#endif
