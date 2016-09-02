#include "StrainLin.hpp"
using namespace Eigen;

StrainLin::StrainLin()
{
	param.resize(2);
	param[0] = 1;
	param[1] = 10;
}

cfgScalar StrainLin::getEnergy(const Matrix3S & F)
{
	Matrix3S I = Matrix3S::Identity();
	Matrix3S eps = 0.5*(F + F.transpose()) - I;
	cfgScalar t = eps.trace();
	cfgScalar Psi = param[0]*eps.squaredNorm() + 0.5f*param[1] * t*t;
	return Psi;
}

Matrix3S StrainLin::getPK1(const Matrix3S & F)
{
	cfgScalar mu = param[0];
	cfgScalar lambda = param[1];
	Matrix3S I = Matrix3S::Identity();
	Matrix3S PP = mu*(F + F.transpose()-2*I) + lambda * (F.trace()-3) * I;
	return PP;
}

Matrix3S StrainLin::getdPdx(const Matrix3S & F, const Matrix3S & dF)
{
	Matrix3S dP = Matrix3S();
	cfgScalar mu = param[0];
	cfgScalar lambda = param[1];
	Matrix3S I = Matrix3S::Identity();
	dP = mu * (dF + dF.transpose()) + lambda * dF.trace() * I;
	return dP;
}

MatrixXS StrainLin::EMatrix()
{
	MatrixXS E=MatrixXS::Zero(6,6);
	cfgScalar G = param[0];
  if (G==0)
    return E;

	cfgScalar l = param[1];
	cfgScalar e = G*(3 * l + 2 * G) / (l + G);
	cfgScalar nu = l / (2 * (l + G));
	cfgScalar c = e / ((1 + nu)*(1 - 2 * nu));
	for (int ii = 0; ii < 3; ii++){
		for (int jj = 0; jj < 3; jj++){
			E(ii, jj) = nu;
		}
		E(ii, ii) = 1 - nu;
		E(3 + ii, 3 + ii) = 0.5f - nu;
	}
	return c*E;
}

Matrix3S StrainLin::getStrainTensor(const Matrix3S & F)
{
  Matrix3S I = Matrix3S::Identity();
	Matrix3S eps = 0.5*(F + F.transpose()) - I;
  return eps;
}
