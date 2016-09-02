#ifndef STRAINENE_HPP
#define STRAINENE_HPP
#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"
///@brief abstract class for strain energy functions, i.e. function of deformation gradient
///also computes derivatives
class StrainEne{
public:
  virtual cfgScalar getEnergy(const Matrix3S & F) = 0;
  virtual Matrix3S getPK1(const Matrix3S & F) = 0;
  virtual Matrix3S getdPdx(const Matrix3S & F, const Matrix3S & dF) = 0;
  virtual ~StrainEne();
  virtual Matrix3S getStrainTensor(const Matrix3S & F) = 0;
  virtual MatrixXS EMatrix()=0;
  
  std::vector<cfgScalar> param;
};
#endif
