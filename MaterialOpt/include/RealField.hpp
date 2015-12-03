#ifndef REAL_FIELD_HPP
#define REAL_FIELD_HPP

#include <Eigen/Sparse>

class RealField{
public:

    RealField(){}

    virtual double f(const Eigen::VectorXd & x)=0;

    ///@brief df/dparam at x.
    virtual Eigen::SparseVector<double> df(const Eigen::VectorXd & x){
      return Eigen::SparseVector<double>(0);
    }
    virtual void setParam(int paramIdx, double val){}
};


#endif
