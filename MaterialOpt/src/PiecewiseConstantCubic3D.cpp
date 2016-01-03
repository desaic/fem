#include "PiecewiseConstantCubic3D.hpp" 
#include "ArrayUtil.hpp"
void
PiecewiseConstantCubic3D::allocate(int nx, int ny, int nz)
{
  gridSize[0] = nx;
  gridSize[1] = ny;
  gridSize[2] = nz;
  int Nparam = nx*ny*nz;
  param.resize(Nparam, 1);
}

Eigen::VectorXi
PiecewiseConstantCubic3D::gridIdx(const Eigen::VectorXd & x)
{
  int dim = 3;
  //cell index in x and y directions.
  Eigen::VectorXi idx(dim);
  Eigen::VectorXd coord = 2*x;
  for (int dd = 0; dd < dim; dd++){
    //reflect other squares to the bottom left square.
    if (coord[dd] > 1){
      coord[dd] = 2 - coord[dd];
    }    
  }

  //bubble sort 3 axis so that x>=y>=z
  double tmp = coord[0];
  if (coord[0] < coord[1]){
    tmp = coord[0];
    coord[0] = coord[1];
    coord[1] = tmp;
  }
  if (coord[0] < coord[2]){
    tmp = coord[0];
    coord[0] = coord[2];
    coord[2] = tmp;
  }
  if (coord[1] < coord[2]){
    tmp = coord[1];
    coord[1] = coord[2];
    coord[2] = tmp;
  }
  for (int dd = 0; dd < dim; dd++){
    int ii = (int)(coord[dd]*gridSize[dd]);
    ii = std::max(ii, 0);
    ii = std::min(ii, gridSize[dd]-1);
    idx[dd] = ii;
  }
  return idx;
}

double
PiecewiseConstantCubic3D::f(const Eigen::VectorXd & x)
{

  assert(x.rows() >= 3);

  Eigen::VectorXi idx = gridIdx(x);
  double val = param[linearIdx(idx, gridSize)];
  return val;
}

Eigen::SparseVector<double>
PiecewiseConstantCubic3D::df(const Eigen::VectorXd & x)
{
  Eigen::VectorXi idx = gridIdx(x);
  Eigen::SparseVector<double> vec(param.size());
  vec.insert(linearIdx(idx, gridSize)) = 1;
  return vec;
}
