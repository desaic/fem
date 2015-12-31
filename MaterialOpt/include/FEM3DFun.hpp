#ifndef FEM3DFUN_HPP
#define FEM3DFUN_HPP

#include "cfgDefs.h"
#include "RealFun.hpp"
#include <vector>

class ElementMesh;
struct PardisoState;
class RealField;

///@brief a real-valued function based on 2D FEM simulation.
class FEM3DFun :public RealFun{
public:

  FEM3DFun();
  ~FEM3DFun();

  ///@brief not freed by ~FEM2DFun destructor.
  ElementMesh * em;
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
  bool constrained;
  ///@brief 6x6 matrix of coarse strain tensor
  ///produced by harmonic displacements.
  ///Each column is a strain.
  Eigen::MatrixXd G;

  ///@brief target strains corresponding to 6 harmonic displacements.
  Eigen::MatrixXd G0;

  ///@brief weight for each column of G.
  ///initialized in init() to all 1s.
  Eigen::MatrixXd wG;
  
  ///@brief 6x6 fine energy matrix. 
  ///\sum_i G^T:C:G
  //Eigen::MatrixXd GTCG;

  double density;

  ///@brief target density fraction.
  double m0;
  ///@brief weight for density term.
  double mw;

  ///@brief external forces
  std::vector< std::vector<double> > externalForce;
  
  ///@brief total force applied on one side
  double forceMagnitude;

  std::vector<int> gridSize;

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
  
  PardisoState * pardisoState;

  MatrixXS K0;
  MatrixXS getKe(int ei);
  void getStiffnessSparse();

  void init(const Eigen::VectorXd & x0);
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
void stretchX(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int>& gridSize, std::vector<double> & externalForce);
void stretchY(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int>& gridSize, std::vector<double> & externalForce);
void stretchZ(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int>& gridSize, std::vector<double> & externalForce);
void shearXY(ElementMesh * em, double ff, const std::vector<int>& gridSize, std::vector<double> & externalForce);
void shearYZ(ElementMesh * em, double ff, const std::vector<int>& gridSize, std::vector<double> & externalForce);
void shearXZ(ElementMesh * em, double ff, const std::vector<int>& gridSize, std::vector<double> & externalForce);

int gridToLinearIdx(int ix, int iy, int iz, const std::vector<int> & gridSize);
#endif
