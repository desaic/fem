#ifndef REAL_FIELD_HPP
#define REAL_FIELD_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>

class RealField{
public:

  RealField(){}

  ///@param x coordinate in space. Usually 2D or 3D.
  virtual double f(const Eigen::VectorXd & x) = 0;

  ///@brief df/dparam at x.
  virtual Eigen::SparseVector<double> df(const Eigen::VectorXd & x){
    return Eigen::SparseVector<double>(0);
  }

  virtual void setParam(int paramIdx, double val){
    assert(paramIdx >= 0 && paramIdx<param.size());
    param[paramIdx] = val;
  }

  Eigen::VectorXd param;
};


#endif
