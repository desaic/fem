#ifndef BILINEAR_FIELD_2D_HPP
#define BILINEAR_FIELD_2D_HPP

#include "RealField.hpp"

///@brief the grid is column major.
class BilinearField2D : public RealField{
public:
  BilinearField2D() :gridSize(2, 0){}

  ///@brief allocate a nx by ny grid. The number of control points is
  /// (nx+1)x(ny+1).
  void allocate(int nx, int ny);

  ///@param x should be in [0,1]
  virtual double f(const Eigen::VectorXd & x);

  virtual Eigen::SparseVector<double> df(const Eigen::VectorXd & x);

  virtual void setParam(int paramIdx, double val){
    assert(paramIdx>=0 && paramIdx<param.size());
    param[paramIdx] = val;
  }

  int vertexIndex(int xi, int yi) const;

  Eigen::VectorXi gridIdx(const Eigen::VectorXd & x);

  ///@brief natural coordinate of x inside a square [0 1].
  Eigen::VectorXd natCoord(const Eigen::VectorXd & x);

  ///@brief values on vertices of a square.
  Eigen::VectorXd squareVertexVal(const Eigen::VectorXi & idx);

  ///@brief vertex indices of a square.
  Eigen::VectorXi squareVertexIdx(const Eigen::VectorXi & idx);

  std::vector<double> param;

  ///@brief gridSize[0] = number of columns, number of cells in x direction.
  std::vector<int> gridSize;
};


#endif
