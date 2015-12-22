#ifndef STRAIN_LIN_HPP
#define STRAIN_LIN_HPP
#include "StrainEne.hpp"
#include <Eigen/Dense>
class StrainLin :public StrainEne{
public:
	StrainLin();
  virtual float getEnergy(const Eigen::Matrix3f & F);
  virtual Eigen::Matrix3f getPK1(const Eigen::Matrix3f & F);
  virtual Eigen::Matrix3f getdPdx(const Eigen::Matrix3f & F, const Eigen::Matrix3f & dF);

	virtual Eigen::MatrixXf EMatrix();

  virtual Eigen::Matrix3f getStrainTensor(const Eigen::Matrix3f & F);
};

#endif
