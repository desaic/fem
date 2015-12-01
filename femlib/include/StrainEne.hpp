#ifndef STRAINENE_HPP
#define STRAINENE_HPP
#include "vecmath.h"
#include <vector>
#include <Eigen\Dense>
///@brief abstract class for strain energy functions, i.e. function of deformation gradient
///also computes derivatives
class StrainEne{
public:
  virtual float getEnergy(const Matrix3f & F)=0;
  virtual Matrix3f getPK1(const Matrix3f & F)=0;
  virtual Matrix3f getdPdx(const Matrix3f & F,const Matrix3f & dF)=0;
  virtual ~StrainEne();
  virtual Matrix3f getStrainTensor(const Matrix3f & F)=0;
  virtual Eigen::MatrixXf EMatrix()=0;
  
  std::vector<float> param;
};
#endif