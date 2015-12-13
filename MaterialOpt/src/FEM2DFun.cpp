#include "FEM2DFun.hpp"
#include "Element2D.h"
#include "ElementMesh2D.h"
#include "RealField.hpp"
#include "SparseLin.hpp"

typedef Eigen::Triplet<cfgScalar> TripletS;
typedef std::vector<std::vector< int> > Grid2D;

///@brief a copy of stiffness matrix assembly function in ElementMesh2D class.
///Difference is that this function scales K for each element using param.
void getStiffnessSparse(ElementMesh2D * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic);

std::vector<int> topVerts(ElementMesh2D * em, const Grid2D & grid);
std::vector<int> botVerts(ElementMesh2D * em, const Grid2D & grid);
std::vector<int> leftVerts(ElementMesh2D * em, const Grid2D & grid);
std::vector<int> rightVerts(ElementMesh2D * em, const Grid2D & grid);

///@brief make stretching force in x direction
void stretchX(ElementMesh2D * em, const Eigen::Vector3d & ff, const Grid2D& grid, std::vector<double> & externalForce);

double measureStretchX(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid);
double measureStretchY(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid);

void FEM2DFun::init(const Eigen::VectorXd & x0){
  bool triangular = true;
  m_I.clear();
  m_J.clear();
  em->stiffnessPattern(m_I, m_J, triangular, m_fixRigid, m_fixRigid, m_periodic);
  sparseInit();
  param = x0;

  int nrow = m_I.size()-1;
  u.resize(externalForce.size());
  for (unsigned int ii = 0; ii < u.size(); ii++){
    u.resize(nrow);
  }
  dfdu.resize(u.size());

  grid.resize(m_nx);
  for (int ii = 0; ii < m_nx; ii++){
    grid[ii].resize(m_ny);
    for (int jj = 0; jj < m_ny; jj++){
      grid[ii][jj] = ii * m_ny + jj;
    }
  }
}

void FEM2DFun::setParam(const Eigen::VectorXd & x0)
{
  bool triangle = true;
  bool constrained = false;
  param = x0;

  //solve linear statics problems
  std::vector<cfgScalar> val;
  getStiffnessSparse(em, param, val, triangle, constrained, m_fixRigid, m_periodic);
  m_val = std::vector<double>(val.begin(), val.end());
  double nrows = m_I.size() - 1;
  for (unsigned int ii = 0; ii < externalForce.size(); ii++){
    std::fill(u[ii].begin(), u[ii].end(), 0);
    sparseSolve(m_I.data(), m_J.data(), m_val.data(), nrows, u[ii].data(), externalForce[ii].data());
  }
}

double FEM2DFun::f()
{
  //Example: measure width and height under a stretching force.
  double dx, dy;
  dx = measureStretchX(em, u[0], grid);
  dy = measureStretchY(em, u[0], grid);
  double val = dxw * (dx - dx0) * (dx - dx0) + dyw * (dy - dy0) * (dy - dy0);
  return val;
}

void FEM2DFun::compute_dfdu()
{
  int nrows = m_I.size() - 1;
  double dx = measureStretchX(em, u[0], grid);
  double dy = measureStretchY(em, u[0], grid);
  std::vector<double>grad(u[0].size(), 0);

  for (unsigned int ii = 0; ii < dfdu.size(); ii++){
    std::vector<int> verts;
    verts = topVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      grad[dim * verts[ii] + 1] += (dy - dy0)*dyw;
    }
    verts = botVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      grad[dim * verts[ii] + 1] -= (dy - dy0)*dyw;
    }
    verts = rightVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      grad[dim * verts[ii]] += (dx - dx0)*dxw;
    }
    verts = leftVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      grad[dim * verts[ii]] -= (dx - dx0)*dxw;
    }
    for (unsigned int ii = 0; ii<grad.size(); ii++){
      grad[ii] /= verts.size();
    }
  }
  dfdu.push_back(grad);
}

Eigen::VectorXd FEM2DFun::df()
{
  Eigen::VectorXd grad(param.size(), 0);
  int nrows = m_I.size() - 1;
  //sensitivity analysis using the adjoint method.
  for (unsigned int ii = 0; ii < u.size(); ii++){
    std::vector<double> lambda(nrows, 0);
    //lambda = K^{-1} dfdu.
    //dfdx = -lambda * dK/dparam * u.
    sparseSolve(m_I.data(), m_J.data(), m_val.data(), nrows, lambda.data(), dfdu[ii].data());
    for (unsigned int jj = 0; jj < param.size(); jj++){
      Element2D * ele = em->e[jj];
      int nV = ele->nV();
      Eigen::MatrixXd dKdp;
      Eigen::VectorXd Ue(nV);
      Eigen::VectorXd lambda_e(nV);
      //in this example dK/dParam = K.
      dKdp = em->getStiffness(jj).cast<double>();
      Ue = dKdp * Ue;
      for (int kk = 0; kk < ele->nV(); kk++){
        int vidx = em->e[jj]->at(kk);
        for (int ll = 0; ll < dim; ll++){
          Ue[kk*dim + ll] = u[ii][vidx*dim + ll];
          lambda_e[kk*dim + ll] = lambda[vidx*dim + ll];
        }
      }
      grad[jj] += -lambda_e.dot(dKdp * Ue);
    }
  }
  return grad;
}

FEM2DFun::FEM2DFun() :em(0), dim(2),
m_Init(false),
m_periodic(true),
m_fixRigid(true),
dx0(1e-3), dy0(1e-3),
dxw(5), dyw(1),
m_nx(0), m_ny(0),
field(0)
{
}

FEM2DFun::~FEM2DFun(){}

void getStiffnessSparse(ElementMesh2D * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic)
{
  int N = 2 * (int)em->x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N, N);
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    Element2D * ele = em->e[ii];
    int nV = ele->nV();
    MatrixXS K = em->getStiffness(ii);

    //scale by parameter.
    //change for more complex material mixture.
    K *= param[ii];

    for (int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for (int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for (int dim1 = 0; dim1<2; dim1++){
          for (int dim2 = 0; dim2<2; dim2++){
            if (trig && (2 * vk + dim2 > 2 * vj + dim1)) {
              continue;
            }
            cfgScalar val = K(2 * jj + dim1, 2 * kk + dim2);
            if (constrained){
              if ( (em->fixed[vk] || em->fixed[vj]) 
                && (vj != vk || dim1 != dim2)){
                val = 0;
              }
            }
            TripletS triple(2 * vj + dim1, 2 * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  if (iFixedRigid)
  {
    em->fixTranslation(Ksparse, trig, em);
    em->fixRotation(Ksparse, trig, em);
  }
  if (iPeriodic)
  {
    em->enforcePeriodicity(Ksparse, trig, em);
  }

  for (int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
      val.push_back(it.value());
    }
  }
}

std::vector<int> topVerts(ElementMesh2D * em, const Grid2D & grid)
{
  std::vector<int> v;
  int nx = grid.size();
  int ny = grid[0].size();
  int topV[2] = { 1, 3 };
  for (int ii = 0; ii<nx; ii++){
    int ei = grid[ii][ny - 1];
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(topV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

std::vector<int> botVerts(ElementMesh2D * em, const Grid2D & grid)
{
  std::vector<int> v;
  int nx = grid.size();
  int botV[2] = { 0, 2 };
  for (int ii = 0; ii<nx; ii++){
    int ei = grid[ii][0];
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(botV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

std::vector<int> leftVerts(ElementMesh2D * em, const Grid2D & grid)
{
  std::vector<int> v;
  int ny = grid[0].size();
  int leftV[2] = { 0, 1 };
  for (int ii = 0; ii<ny; ii++){
    int ei = grid[0][ii];
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(leftV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

std::vector<int> rightVerts(ElementMesh2D * em, const Grid2D & grid)
{
  std::vector<int> v;
  int nx = grid.size();
  int ny = grid[0].size();
  int rightV[2] = { 2, 3 };
  for (int ii = 0; ii<ny; ii++){
    int ei = grid[nx - 1][ii];
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(rightV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

void stretchX(ElementMesh2D * em, const Eigen::Vector3d & ff, const Grid2D& grid, std::vector<double> & externalForce)
{
  int dim = 2;
  std::vector<int> leftv, rightv;
  rightv = rightVerts(em, grid);
  Eigen::Vector3d fv = ff / rightv.size();
  for (unsigned int ii = 0; ii<rightv.size(); ii++){
    int vidx = rightv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx+jj] += fv[jj];
    }
  }
  leftv = leftVerts(em, grid);
  for (unsigned int ii = 0; ii<leftv.size(); ii++){
    int vidx = leftv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx + jj] -= fv[jj];
    }
  }
}

double measureStretchX(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid)
{
  int dim = 2;
  std::vector<int> lv, rv;
  lv = leftVerts(em, grid);
  rv = rightVerts(em, grid);
  double stretch = 0;
  for (unsigned int ii = 0; ii<lv.size(); ii++){
    stretch += u[dim * rv[ii]] - u[dim * lv[ii]];
  }
  stretch /= lv.size();
  return stretch;
}

double measureStretchY(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid)
{
  int dim = 2;
  std::vector<int> tv, bv;
  tv = topVerts(em, grid);
  bv = botVerts(em, grid);
  double stretch = 0;
  for (unsigned int ii = 0; ii<tv.size(); ii++){
    stretch += u[dim * tv[ii]+1] - u[dim * bv[ii]+1];
  }
  stretch /= bv.size();
  return stretch;
}
