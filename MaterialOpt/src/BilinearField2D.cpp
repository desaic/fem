#include "BilinearField2D.hpp" 

//1 3
//0 2
Eigen::VectorXd bilinearWeight(const Eigen::VectorXd & n){
  Eigen::VectorXd val(4);
  val[0] = (1-n[0]) * (1-n[1]);
  val[1] = (1-n[0]) * (  n[1]);
  val[2] = (  n[0]) * (1-n[1]);
  val[3] = (  n[0]) * (  n[1]);
  return val;
}

void
BilinearField2D::allocate(int nx, int ny)
{
  gridSize[0] = nx;
  gridSize[1] = ny;
  int Nparam = (nx+1)*(ny+1);
  param.resize(Nparam, 0);
}

Eigen::VectorXi
BilinearField2D::gridIdx(const Eigen::VectorXd & x)
{
  int dim = 2;
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

Eigen::VectorXd
BilinearField2D::natCoord(const Eigen::VectorXd & x)
{
  int dim = 2;
  Eigen::VectorXd n(dim);
  Eigen::VectorXi idx = gridIdx(x);
  for(int dd = 0; dd<dim; dd++){
    double val = x[dd] * gridSize[dd] - idx[dd];
    val = std::max(0.0, val);
    val = std::min(1.0, val);
    n[dd] = val;
  }
  return n;
}

Eigen::VectorXi
BilinearField2D::squareVertexIdx(const Eigen::VectorXi & idx)
{
  int nV = 4;
  Eigen::VectorXi val(nV);
  int vi[4][2] ={{0,0},{0,1},{1,0},{1,1}};
  for(int ii = 0; ii<nV; ii++){
    val[ii] = vertexIndex(idx[0] + vi[ii][0], idx[1] + vi[ii][1]);
  }
  return val;
}

Eigen::VectorXd
BilinearField2D::squareVertexVal(const Eigen::VectorXi & idx)
{
  Eigen::VectorXi vi = squareVertexIdx(idx);
  int nV = vi.rows();
  Eigen::VectorXd val(nV);
  for(int ii = 0; ii<nV; ii++){
    val[ii] = param[vi[ii]];
  }
  return val;
}

int
BilinearField2D::vertexIndex(int xi, int yi) const
{
  return xi*(gridSize[1]+1) + yi;
}

double
BilinearField2D::f(const Eigen::VectorXd & x)
{

  assert(x.rows() >= 2);

  Eigen::VectorXd n = natCoord(x);
  Eigen::VectorXd w = bilinearWeight(n);
  Eigen::VectorXi idx = gridIdx(x);
  Eigen::VectorXd val= squareVertexVal(idx);
  return val.dot(w);
}

Eigen::SparseVector<double>
BilinearField2D::df(const Eigen::VectorXd & x)
{

  Eigen::VectorXi idx = gridIdx(x);
  Eigen::VectorXi vi = squareVertexIdx(idx);
  Eigen::VectorXd n = natCoord(x);
  Eigen::VectorXd w = bilinearWeight(n);
  int nV = vi.rows();
  Eigen::SparseVector<double> vec(param.size());
  vec.reserve(nV);
  for(int ii = 0; ii<nV; ii++){
    vec.insert(vi[ii]) = w[ii];
  }
  return vec;
}
