#ifndef STRAINENE_HPP
#define STRAINENE_HPP
#include <vector>
#include <Eigen/Dense>
///@brief abstract class for strain energy functions, i.e. function of deformation gradient
///also computes derivatives
class StrainEne{
public:
  virtual float getEnergy(const Eigen::Matrix3f & F) = 0;
  virtual Eigen::Matrix3f getPK1(const Eigen::Matrix3f & F) = 0;
  virtual Eigen::Matrix3f getdPdx(const Eigen::Matrix3f & F, const Eigen::Matrix3f & dF) = 0;
  virtual ~StrainEne();
  virtual Eigen::Matrix3f getStrainTensor(const Eigen::Matrix3f & F) = 0;
  virtual Eigen::MatrixXf EMatrix()=0;
  
  std::vector<float> param;
};
#endif
