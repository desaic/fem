#include "StrainEneNeo2D.h"
using namespace Eigen;

StrainEneNeo2D::StrainEneNeo2D()
{
  param.resize(2);
  param[0] =1;
  param[1] =10;
}

cfgScalar StrainEneNeo2D::getEnergy(const Matrix2S & F)
{
  cfgScalar lambda3 = 1;
  cfgScalar traceC = (F.transpose()*F).trace() + lambda3;
  cfgScalar detF = F.determinant();
  cfgScalar I1 = traceC;
  cfgScalar JJ = std::log(detF);
  cfgScalar mu = param[0],lambda=param[1];
  cfgScalar Psi = (mu/2) * (I1-3) - mu*JJ + (lambda/2)*JJ*JJ;
  return (cfgScalar)Psi;
}

Matrix2S StrainEneNeo2D::getPK1(const Matrix2S & F) 
{
  cfgScalar JJ = std::log(F.determinant());
  Matrix2S Finv = F.inverse().transpose();
  cfgScalar mu = param[0],lambda=param[1];
  Matrix2S PP = mu*(F-Finv) + lambda*JJ*Finv;
  return PP;
}

Matrix2S StrainEneNeo2D::getdPdx(const Matrix2S & F,const Matrix2S & dF)
{
  Matrix2S dP = Matrix2S();
  cfgScalar JJ = std::log(F.determinant());
  Matrix2S Finv = F.inverse();
  Matrix2S FinvT = Finv.transpose();
  cfgScalar mu = param[0],lambda=param[1];
  dP = mu*dF;
  cfgScalar c1 = mu-lambda * JJ;
  dP += c1 * FinvT*dF.transpose()*FinvT;
  dP += lambda*(Finv*dF).trace()*FinvT;
  return dP;
}

Matrix2S StrainEneNeo2D::getStressTensor(const Matrix2S & F)
{
  assert(0); // not implemented
  return Matrix2S::Identity();
}

Matrix2S StrainEneNeo2D::getStrainTensor(const Matrix2S & F)
{
  assert(0); // not implemented
  return Matrix2S::Identity();
}

MatrixXS StrainEneNeo2D::getElasticityTensor()
{
  assert(0); // not implemented
  return Matrix2S::Identity();
}

