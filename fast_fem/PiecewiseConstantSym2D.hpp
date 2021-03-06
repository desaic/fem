#ifndef PIECEWISE_CONSTANT_SYM_2D_HPP
#define PIECEWISE_CONSTANT_SYM_2D_HPP

#include "RealField.hpp"

///@brief the grid is column major.
///Symmetric version of PiecewiseConstant2D
class PiecewiseConstantSym2D : public RealField{
public:
  PiecewiseConstantSym2D() :gridSize(2, 0){}

  ///@brief allocate a nx by ny grid. The number of control points is
  /// (nx+1)x(ny+1).
  void allocate(int nx, int ny);

  ///@param x should be in [0,1]
  virtual double f(const Eigen::VectorXd & x);

  virtual Eigen::SparseVector<double> df(const Eigen::VectorXd & x);

  int vertexIndex(int xi, int yi) const;

  ///@brief grid index given a point x.
  ///@param x should be inside [0,1]^2.
  ///This is the only function different from PiecewiseConstant2D
  Eigen::VectorXi gridIdx(const Eigen::VectorXd & x);

  ///@brief gridSize[0] = number of columns, number of cells in x direction.
  std::vector<int> gridSize;
};


#endif
