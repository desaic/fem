#ifndef STRAINCOROTLIN_HPP
#define STRAINCOROTLIN_HPP
#include "StrainEne.hpp"
class StrainCorotLin :public StrainEne{
public:
  StrainCorotLin();
  virtual float getEnergy(const Matrix3f & F);
  virtual Matrix3f getPK1(const Matrix3f & F);
  virtual Matrix3f getdPdx(const Matrix3f & F, const Matrix3f & dF);
  
};

#endif