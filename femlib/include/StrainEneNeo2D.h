#ifndef StrainEneNeo2D_h
#define StrainEneNeo2D_h
#include "StrainEne2D.h"

class StrainEneNeo2D:public StrainEne2D
{
public:
  StrainEneNeo2D();
  virtual cfgScalar getEnergy(const Matrix2S & F);
  virtual Matrix2S getPK1(const Matrix2S & F);
  virtual Matrix2S getdPdx(const Matrix2S & F, const Matrix2S & dF);

  virtual Matrix2S getStressTensor(const Matrix2S & F);
  virtual Matrix2S getStrainTensor(const Matrix2S & F);
  virtual MatrixXS getElasticityTensor();

};
#endif