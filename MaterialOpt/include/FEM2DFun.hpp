#ifndef FEM2DFUN_HPP
#define FEM2DFUN_HPP

#include "RealFun.hpp"
#include <vector>

class ElementMesh2D;
class RealField;

///@brief a real-valued function based on 2D FEM simulation.
class FEM2DFun :public RealFun{
public:

  FEM2DFun();
  ~FEM2DFun();

  ///@brief not freed by ~FEM2DFun destructor.
  ElementMesh2D * em;
  ///@brief dimension of space. set to 2.
  int dim;
  ///@brief variables used by simulation and linear solve.
  ///@brief indices of stiffness matrix.
  ///@brief m_I is array of length matrix_rows + 1.
  std::vector<int> m_I, m_J;
  ///@brief values in stiffness matrix. Updated by setParam call.
  std::vector<double> m_val;
  bool m_Init;
  bool m_periodic;
  bool m_fixRigid;

  ///@brief target displacements.
  double dx0, dy0;
  ///@brief weight for displacement objectives.
  double dxw, dyw;

  ///@brief external forces
  std::vector< std::vector<double> > externalForce;
  ///@brief total force applied on one side
  double forceMagnitude;

  ///@brief element indices ordered into a grid
  ///grid[col][row] is the element index.
  std::vector < std::vector< int> > grid;
  int m_nx, m_ny;

  ///@brief linear static displacements solved using external forces.
  std::vector< std::vector<double> > u;
  
  ///@brief gradient of objective with respect to displacements.
  ///Each vector should be the same size as the stiffness matrix, padded by 0s
  ///so that the call to linear solver works.
  std::vector< std::vector<double> > dfdu;

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

  ///@brief this function needs to change according to f().
  void compute_dfdu();

  ///@brief compute gradient of f().
  ///Given an implementation of dfdu.
  ///When dK/dParam changes, this function needs to change accordingly.
  ///The return value should have the same size as the parameters.
  virtual Eigen::VectorXd df();

};

#endif
