#include "StrainLin2D.h"

StrainLin2D::StrainLin2D()
{
	param.resize(2);
	param[0] = 1;
	param[1] = 10;
}

cfgScalar StrainLin2D::getEnergy(const Matrix2S & F)
{
	Matrix2S I = Matrix2S::Identity();
	Matrix2S eps = 0.5*(F + F.transpose()) - I;
	cfgScalar t = eps.trace();
  cfgScalar Psi = param[0]*(eps(0,0)*eps(0,0) + eps(0,1)*eps(0,1) + eps(1,0)*eps(1,0) + eps(1,1)*eps(1,1)) + 0.5f*param[1] * t*t;
	//cfgScalar Psi = param[0]*eps.norm2() + 0.5f*param[1] * t*t;
	return Psi;
}

Matrix2S StrainLin2D::getPK1(const Matrix2S & F)
{
	cfgScalar mu = param[0];
	cfgScalar lambda = param[1];
	Matrix2S I = Matrix2S::Identity();
	Matrix2S PP = mu*(F + F.transpose()-2*I) + lambda * (F.trace()-2) * I;
	return PP;
}

Matrix2S StrainLin2D::getStrainTensor(const Matrix2S & F)
{
  Matrix2S I = Matrix2S::Identity();
	Matrix2S eps = 0.5*(F + F.transpose()) - I;
  return eps;
}

Matrix2S StrainLin2D::getStressTensor(const Matrix2S & F)
{
  cfgScalar & mu = param[0];
  cfgScalar & lambda = param[1];

  Matrix2S I = Matrix2S::Identity();
	Matrix2S eps = 0.5*(F + F.transpose()) - I;
  
  Matrix2S s = 2*mu*eps + lambda*eps.trace()*I;
  return s;
}

MatrixXS StrainLin2D::getElasticityTensor()
{
  MatrixXS C = MatrixXS::Zero(3, 3);

  cfgScalar & mu = param[0];
  cfgScalar & lambda = param[1];

  cfgScalar E = mu*(3*lambda+2*mu)/(lambda+mu);
  cfgScalar nu = lambda/(2*(lambda+mu));

  // plane stress
  cfgScalar coeff = E/(1-nu*nu);
  C(0,0) = 1;
  C(0,1) = nu;
  C(1,0) = nu;
  C(1,1) = 1;
  C(2,2) = (1-nu)/2;
  C *= coeff;

  // plane strain
  /*C(0,0) = 2*mu+lambda;
  C(0,1) = lambda;
  C(1,0) = lambda;
  C(1,1) = 2*mu+lambda;
  C(2,2) = mu;*/ 

  return C;
}

Matrix2S StrainLin2D::getdPdx(const Matrix2S & F, const Matrix2S & dF)
{
	Matrix2S dP = Matrix2S();
	cfgScalar mu = param[0];
	cfgScalar lambda = param[1];
	Matrix2S I = Matrix2S::Identity();
	dP = mu * (dF + dF.transpose()) + lambda * dF.trace() * I;
	return dP;
}
