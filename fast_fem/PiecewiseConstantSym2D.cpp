#include "PiecewiseConstantSym2D.hpp" 

void
PiecewiseConstantSym2D::allocate(int nx, int ny)
{
  gridSize[0] = nx;
  gridSize[1] = ny;
  int Nparam = nx*ny;
  param.resize(Nparam, 1);
}

Eigen::VectorXi
PiecewiseConstantSym2D::gridIdx(const Eigen::VectorXd & x)
{
  int dim = 2;
  //cell index in x and y directions.
  Eigen::VectorXi idx(dim);
  for(int dd = 0; dd<dim; dd++){
    double coord = 2*x[dd];
    //reflect other squares to the bottom left square.
    if (coord > 1){
      coord = 2 - coord;
    }
    int ii = (int)(coord*gridSize[dd]);
    ii = std::max(ii, 0);
    ii = std::min(ii, gridSize[dd]-1);
    idx[dd] = ii;
  }
  return idx;
}

double
PiecewiseConstantSym2D::f(const Eigen::VectorXd & x)
{

  assert(x.rows() >= 2);

  Eigen::VectorXi idx = gridIdx(x);
  double val= param[idx[0] * gridSize[1] + idx[1]];
  return val;
}

Eigen::SparseVector<double>
PiecewiseConstantSym2D::df(const Eigen::VectorXd & x)
{
  Eigen::VectorXi idx = gridIdx(x);
  Eigen::SparseVector<double> vec(param.size());
  int vidx = idx[0] * gridSize[1] + idx[1];
  vec.insert(vidx) = 1;
  return vec;
}
