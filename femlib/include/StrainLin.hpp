#ifndef STRAIN_LIN_HPP
#define STRAIN_LIN_HPP
#include "StrainEne.hpp"
#include <vecmath.h>
#include <Eigen\Dense>
class StrainLin :public StrainEne{
public:
	StrainLin();
	float getEnergy(const Matrix3f & F);
	Matrix3f getPK1(const Matrix3f & F);
	Matrix3f getdPdx(const Matrix3f & F, const Matrix3f & dF);

	Eigen::MatrixXf EMatrix();
};

#endif