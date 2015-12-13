#ifndef FEM2DFUN_HPP
#define FEM2DFUN_HPP

#include "RealFun.hpp"

class ElementMesh2D;
class RealField;

///@brief a real-valued function based on 2D FEM simulation.
class FEM2DFun :public RealFun{
public:

  FEM2DFun();
  ~FEM2DFun();

  ///@brief not freed by ~FEM2DFun destructor.
  ElementMesh2D * em;

  ///@brief variables used by simulation and linear solve.
  std::vector<int> m_I, m_J;
  bool m_Init;
  bool m_periodic;
  bool m_noRigid;

  ///@brief not freed by ~FEM2DFun destructor.
  ///parameterization for material assignments.
  ///BilinearField2D for bilinearly interpolated material distributsions
  ///over a grid and PiecewiseConstant2D for constant material in each cell.
  RealField * field;

  void init(const Eigen::VectorXd & x0);
  ///@brief update parameters for computing function value
  ///and gradients.
  ///In case of fem, it runs required fem simulations.
  ///@param x0 material distributions. Right now assumed to just be one
  ///number for each cell for ratio between two materials.
  virtual void setParam(const Eigen::VectorXd & x0);

  ///@brief evaluate function at the currently set parameter.
  ///Change to optimize for different objective.
  ///The example computes change in width and height under a single
  ///stretching force along x-axis.
  virtual double f();

  ///@brief compute gradient of f().
  ///The return value should have the same size as the parameters.
  virtual Eigen::VectorXd df();

};

#endif
