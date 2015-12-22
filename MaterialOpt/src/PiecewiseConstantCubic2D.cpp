#include "PiecewiseConstantCubic2D.hpp" 

void
PiecewiseConstantCubic2D::allocate(int nx, int ny)
{
  gridSize[0] = nx;
  gridSize[1] = ny;
  int Nparam = nx*ny;
  param.resize(Nparam, 1);
}

Eigen::VectorXi
PiecewiseConstantCubic2D::gridIdx(const Eigen::VectorXd & x)
{
  int dim = 2;
  //cell index in x and y directions.
  Eigen::VectorXi idx(dim);
  Eigen::VectorXd coord = 2*x;
  for (int dd = 0; dd < dim; dd++){
    //reflect other squares to the bottom left square.
    if (coord[dd] > 1){
      coord[dd] = 2 - coord[dd];
    }    
  }
  //reflect along diagonal
  if (coord[0] < coord[1]){
    double tmp = coord[0];
    coord[0] = coord[1];
    coord[1] = tmp;
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
PiecewiseConstantCubic2D::f(const Eigen::VectorXd & x)
{

  assert(x.rows() >= 2);

  Eigen::VectorXi idx = gridIdx(x);
  double val= param[idx[0] * gridSize[1] + idx[1]];
  return val;
}

Eigen::SparseVector<double>
PiecewiseConstantCubic2D::df(const Eigen::VectorXd & x)
{
  Eigen::VectorXi idx = gridIdx(x);
  Eigen::SparseVector<double> vec(param.size());
  int vidx = idx[0] * gridSize[1] + idx[1];
  vec.insert(vidx) = 1;
  return vec;
}
