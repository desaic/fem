#ifndef STRAINENE2D_HPP
#define STRAINENE2D_HPP
#include "vecmath.h"
#include <vector>
///@brief abstract class for strain energy functions, i.e. function of deformation gradient
///also computes derivatives
class StrainEne2D{
public:
  virtual float getEnergy(const Matrix2f & F)=0;
  virtual Matrix2f getPK1(const Matrix2f & F)=0;
  virtual Matrix2f getdPdx(const Matrix2f & F,const Matrix2f & dF)=0;
  virtual ~StrainEne2D() {};
  
  std::vector<float> param;
};
#endif