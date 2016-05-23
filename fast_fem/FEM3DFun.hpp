#ifndef FEM3DFUN_HPP
#define FEM3DFUN_HPP

#include "cfgDefs.h"
#include "RealFun.hpp"
#include "ElasticHexFEM.h"
#include "Homogenization.h"
#include <vector>

class RealField;

///@brief a real-valued function based on 2D FEM simulation.
class FEM3DFun :public RealFun{
public:

  FEM3DFun();
  ~FEM3DFun();
  ///@brief the grid size for lattice is 1 less than gridSize
  ///for periodic boundary conditions.
  Lattice<3> lattice;
  ElasticHexFEM<3> fem;
  Homogenization<3> * hom;
  //displacements. size 6.
  std::vector<Eigen::VectorXd> u;
  //stiffness of base material.
  std::vector<Eigen::MatrixXd> Ke0;

  bool m_Init;
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

  std::vector<int> gridSize;
  
  ///@brief not freed by ~FEM2DFun destructor.
  ///parameterization for material assignments.
  ///BilinearField2D for bilinearly interpolated material distributsions
  ///over a grid and PiecewiseConstant2D for constant material in each cell.
  RealField * field;

  ///@brief material distribution. Currently just ratio between two materials for each element.
  ///Should be the same size as number of elements
  Eigen::VectorXd distribution;
  
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


#endif
