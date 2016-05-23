#ifndef STRAIN_LIN_HPP
#define STRAIN_LIN_HPP
#include "StrainEne.hpp"
#include <Eigen/Dense>
class StrainLin :public StrainEne{
public:
	StrainLin();
  virtual cfgScalar getEnergy(const Matrix3S & F);
  virtual Matrix3S getPK1(const Matrix3S & F);
  virtual Matrix3S getdPdx(const Matrix3S & F, const Matrix3S & dF);

	virtual MatrixXS EMatrix();

  virtual Matrix3S getStrainTensor(const Matrix3S & F);
};

#endif
