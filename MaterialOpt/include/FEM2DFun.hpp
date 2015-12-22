#ifndef FEM2DFUN_HPP
#define FEM2DFUN_HPP

#include "RealFun.hpp"
#include <vector>

class ElementMesh2D;
class RealField;

typedef std::vector<std::vector< int> > Grid2D;

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

  ///@brief displacements produced by the first externalForce after calling setParam(x).
  double dx, dy;
  double density;
  ///@brief target displacements.
  double dx0, dy0;
  ///@brief weight for displacement objectives.
  double dxw, dyw;

  ///@brief target density fraction.
  double m0;
  ///@brief weight for density term.
  double mw;

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

  ///@brief material distribution. Currently just ratio between two materials for each element.
  ///Should be the same size as number of elements
  Eigen::VectorXd distribution;

  void init(const Eigen::VectorXd & x0);
  void computeGrid();
  void initArrays();
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

  ///@brief How element stiffness changes with respect to parameter.
  ///@TODO: implement this so that df() doesn't need to change when dK/dparam changes.
  Eigen::MatrixXd dKdp(int eidx, int pidx);

  ///@brief compute gradient of f().
  ///Given an implementation of dfdu and dKdp.
  ///When dK/dParam changes, this function needs to change accordingly.
  ///The return value should have the same size as the parameters.
  virtual Eigen::VectorXd df();

  void log(std::ostream & out);
};

///@brief make stretching force in x direction
void stretchX(ElementMesh2D * em, const Eigen::Vector3d & ff, const Grid2D& grid, std::vector<double> & externalForce);
void stretchY(ElementMesh2D * em, const Eigen::Vector3d & ff, const Grid2D& grid, std::vector<double> & externalForce);

double measureStretchX(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid);
double measureStretchY(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid);

//measure shear displacement in x direction of top and bottom vertices.
double measureShearX(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid);

#endif
