#ifndef STRAINENENEO_HPP
#define STRAINENENEO_HPP
#include "StrainEne.hpp"

class StrainEneNeo:public StrainEne
{
public:
  StrainEneNeo();
  float getEnergy(const Eigen::Matrix3f & F);
  Eigen::Matrix3f getPK1(const Eigen::Matrix3f & F);
  Eigen::Matrix3f getdPdx(const Eigen::Matrix3f & F, const Eigen::Matrix3f & dF);

  Eigen::MatrixXf EMatrix();
  Eigen::Matrix3f getStrainTensor(const Eigen::Matrix3f & F);
};
#endif