#ifndef STRAINENE2D_HPP
#define STRAINENE2D_HPP
#include "vecmath.h"
#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"

///@brief abstract class for strain energy functions, i.e. function of deformation gradient
///also computes derivatives
class StrainEne2D{
public:
  virtual cfgScalar getEnergy(const Matrix2S & F)=0;
  virtual Matrix2S getPK1(const Matrix2S & F)=0;
  virtual Matrix2S getdPdx(const Matrix2S & F,const Matrix2S & dF)=0;
  virtual Matrix2S getStressTensor(const Matrix2S & F) = 0;
  virtual Matrix2S getStrainTensor(const Matrix2S & F) = 0;
  virtual MatrixXS getElasticityTensor() = 0;

  virtual ~StrainEne2D() {};
  
  std::vector<cfgScalar> param;
};
#endif