#ifndef STRAIN_LIN_2D_HPP
#define STRAIN_LIN_2D_HPP
#include "StrainEne2D.h"
#include <vecmath.h>
#include <Eigen\Dense>
class StrainLin2D :public StrainEne2D{
public:
	StrainLin2D();
  virtual float getEnergy(const Matrix2f & F);
  virtual Matrix2f getPK1(const Matrix2f & F);
  virtual Matrix2f getdPdx(const Matrix2f & F, const Matrix2f & dF);

	virtual Eigen::MatrixXf EMatrix();
};

#endif