#include "PiecewiseConstant3D.hpp" 

void
PiecewiseConstant3D::allocate(int nx, int ny, int nz)
{
  gridSize[0] = nx;
  gridSize[1] = ny;
  gridSize[2] = nz;
  int Nparam = nx*ny*nz;
  param.resize(Nparam, 1);
}

Eigen::VectorXi
PiecewiseConstant3D::gridIdx(const Eigen::VectorXd & x)
{
  int dim = 3;
  //cell index in x and y directions.
  Eigen::VectorXi idx(dim);
  for(int dd = 0; dd<dim; dd++){
    int ii = (int)(x[dd]*gridSize[dd]);
    ii = std::max(ii, 0);
    ii = std::min(ii, gridSize[dd]-1);
    idx[dd] = ii;
  }
  return idx;
}

int linearIdx(const Eigen::VectorXi & idx, const std::vector<int> & gridSize)
{
	return idx[0] * gridSize[1] * gridSize[2] + idx[1] * gridSize[2] + idx[2];
}

double
PiecewiseConstant3D::f(const Eigen::VectorXd & x)
{

  assert(x.rows() >= 3);

  Eigen::VectorXi idx = gridIdx(x);
  double val= param[linearIdx(idx, gridSize)];
  return val;
}

Eigen::SparseVector<double>
PiecewiseConstant3D::df(const Eigen::VectorXd & x)
{
  Eigen::VectorXi idx = gridIdx(x);
  Eigen::SparseVector<double> vec(param.size());
  vec.insert(linearIdx(idx, gridSize)) = 1;
  return vec;
}
