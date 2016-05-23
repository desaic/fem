#ifndef STRAINENENEO_HPP
#define STRAINENENEO_HPP
#include "StrainEne.hpp"

class StrainEneNeo:public StrainEne
{
public:
  StrainEneNeo();
  cfgScalar getEnergy(const Matrix3S & F);
  Matrix3S getPK1(const Matrix3S & F);
  Matrix3S getdPdx(const Matrix3S & F, const Matrix3S & dF);

  MatrixXS EMatrix();
  Matrix3S getStrainTensor(const Matrix3S & F);
};
#endif