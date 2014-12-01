#include "StrainLin.hpp"

StrainLin::StrainLin()
{
	param.resize(2);
	param[0] = 1;
	param[1] = 10;
}

float StrainLin::getEnergy(const Matrix3f & F)
{
	Matrix3f I = Matrix3f::identity();
	Matrix3f eps = 0.5*(F + F.transposed()) - I;
	float t = eps.trace();
	float Psi = param[0]*eps.norm2() + 0.5*param[1] * t*t;
	return (float)Psi;
}

Matrix3f StrainLin::getPK1(const Matrix3f & F)
{
	float mu = param[0];
	float lambda = param[1];
	Matrix3f I = Matrix3f::identity();
	Matrix3f PP = mu*(F + F.transposed()-2*I) + lambda * (F.trace()-3) * I;
	return PP;
}

Matrix3f StrainLin::getdPdx(const Matrix3f & F, const Matrix3f & dF)
{
	Matrix3f dP = Matrix3f();
	float mu = param[0];
	float lambda = param[1];
	Matrix3f I = Matrix3f::identity();
	dP = mu * (dF + dF.transposed()) + lambda * dF.trace() * I;
	return dP;
}
