#ifndef STRAIN_LIN_2D_HPP
#define STRAIN_LIN_2D_HPP
#include "StrainEne2D.h"
#include <vecmath.h>
#include <Eigen\Dense>
class StrainLin2D :public StrainEne2D{
public:
	StrainLin2D();
  virtual cfgScalar getEnergy(const Matrix2S & F);
  virtual Matrix2S getPK1(const Matrix2S & F);
  virtual Matrix2S getdPdx(const Matrix2S & F, const Matrix2S & dF);
  virtual Matrix2S getStressTensor(const Matrix2S & F);
  virtual Matrix2S getStrainTensor(const Matrix2S & F);
  virtual MatrixXS getElasticityTensor();
};

#endif