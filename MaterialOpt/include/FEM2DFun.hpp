#ifndef FEM2DFUN_HPP
#define FEM2DFUN_HPP

#include "RealFun.hpp"

class ElementMesh2D;

class FEM2DFun :public RealFun{
public:
  
  ElementMesh2D * em;

  void init(const Eigen::VectorXd & x0){
    param = x0;
  }

  ///@brief update parameters for computing function value
  ///and gradients. Override to update additional data structures.
  virtual void setParam(const Eigen::VectorXd & x0);

  ///@brief evaluate function at the currently set parameter.
  virtual double f() = 0;

  ///@brief compute gradient. Not implemented by default.
  ///The return value should have the same dimensions as the parameters.
  virtual Eigen::VectorXd df(){
    std::cout << "Unimplemented.\n";
    return Eigen::VectorXd(1);
  }

};


#endif
