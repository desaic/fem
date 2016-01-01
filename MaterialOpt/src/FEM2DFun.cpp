#include "FEM2DFun.hpp"
#include "Element2D.h"
#include "ElementMesh2D.h"
#include "RealField.hpp"
#include "SparseLin.hpp"
//#include "Timer.hpp"

typedef Eigen::Triplet<cfgScalar> TripletS;

///@brief a copy of stiffness matrix assembly function in ElementMesh2D class.
///Difference is that this function scales K for each element using param.
void getStiffnessSparse(ElementMesh2D * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic);

std::vector<int> topVerts(ElementMesh2D * em, const Grid2D & grid);
std::vector<int> botVerts(ElementMesh2D * em, const Grid2D & grid);
std::vector<int> leftVerts(ElementMesh2D * em, const Grid2D & grid);
std::vector<int> rightVerts(ElementMesh2D * em, const Grid2D & grid);

void FEM2DFun::computeGrid()
{
  grid.resize(m_nx);
  for (int ii = 0; ii < m_nx; ii++){
    grid[ii].resize(m_ny);
    for (int jj = 0; jj < m_ny; jj++){
      grid[ii][jj] = ii * m_ny + jj;
    }
  }
}

void FEM2DFun::initArrays()
{
  int nrow = externalForce[0].size();
  u.resize(externalForce.size());
  dfdu.resize(u.size());
  for (unsigned int ii = 0; ii < u.size(); ii++){
    u[ii].resize(nrow);
    dfdu[ii].resize(nrow);
  }
}

void FEM2DFun::init(const Eigen::VectorXd & x0)
{  
  bool triangular = true;
  m_I.clear();
  m_J.clear();
  em->stiffnessPattern(m_I, m_J, triangular, m_fixRigid, m_fixRigid, m_periodic);
  sparseInit();
  param = x0;
  distribution = Eigen::VectorXd::Zero(em->e.size());
  int nrow = (int)m_I.size() - 1;
  if (grid.size() == 0){
    computeGrid();
  }
  externalForce.resize(1);
  externalForce[0].resize(nrow, 0);
  Eigen::Vector3d forceDir = forceMagnitude * Eigen::Vector3d(1, 0, 0);
  stretchX(em, forceDir, grid, externalForce[0]);
  initArrays();
  m_Init = true;
}

void FEM2DFun::setParam(const Eigen::VectorXd & x0)
{
  bool triangle = true;
  bool constrained = false;
  Vector3S color0(0.6, 0.6, 1.0);
  param = x0;

  //read material distribution from field
  for (int ii = 0; ii < field->param.size(); ii++){
    field->setParam(ii, x0[ii]);
  }
  for (int ii = 0; ii < m_nx; ii++){
    for (int jj = 0; jj < m_ny; jj++){
      Eigen::VectorXd coord = Eigen::VectorXd::Zero(dim);
      int eIdx = grid[ii][jj];
      coord[0] = (ii + 0.5) / m_nx;
      coord[1] = (jj + 0.5) / m_ny;
      distribution[eIdx] = field->f(coord);
      em->e[eIdx]->color = distribution[eIdx] * color0;
    }
  }

  //solve linear statics problems
  //Timer timer;
  std::vector<cfgScalar> val;
  //timer.startWall();
  getStiffnessSparse(em, distribution, val, triangle, constrained, m_fixRigid, m_periodic);
  //timer.endWall();
  //std::cout << "assemble time " << timer.getSecondsWall() << "\n";
  m_val = std::vector<double>(val.begin(), val.end());
  int nrows = (int)m_I.size() - 1;
  //checkSparseIndex(m_I, m_J);
  for (unsigned int ii = 0; ii < externalForce.size(); ii++){
    std::fill(u[ii].begin(), u[ii].end(), 0);
    std::vector<int> I = m_I;
    std::vector<int> J = m_J;
    //timer.startWall();
    sparseSolve(I.data(), J.data(), m_val.data(), nrows, &(u[ii][0]), externalForce[ii].data());
    //timer.endWall();
    //std::cout<<"Lin solve time " << timer.getSecondsWall() << "\n";
  }

  dx = measureStretchX(em, u[0], grid);
  dy = measureStretchY(em, u[0], grid);

  //show rendering
  /*for (unsigned int ii = 0; ii < em->x.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      em->x[ii][jj] = em->X[ii][jj] + u[0][ii*dim + jj];
    }
  }*/ 

  //density objective
  density = 0;
  for (unsigned int ii = 0; ii < distribution.size(); ii++){
    density += distribution[ii];
  }
  density /= distribution.size();

}

double FEM2DFun::f()
{
  //Example: measure width and height under a stretching force.
  double val = 0.5 * dxw * (dx - dx0) * (dx - dx0) + 0.5 * dyw * (dy - dy0) * (dy - dy0);
  val += 0.5 * mw * (density - m0) * (density - m0);
  return val ;
}

void FEM2DFun::compute_dfdu()
{
  int nrows = (int)m_I.size() - 1;
  std::fill(dfdu[0].begin(), dfdu[0].end(), 0);
  for (unsigned int ii = 0; ii < dfdu.size(); ii++){
    std::vector<int> verts;
    verts = topVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii] + 1] += (dy - dy0)*dyw;
    }
    verts = botVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii] + 1] -= (dy - dy0)*dyw;
    }
    verts = rightVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii]] += (dx - dx0)*dxw;
    }
    verts = leftVerts(em, grid);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii]] -= (dx - dx0)*dxw;
    }
    for (unsigned int ii = 0; ii<dfdu[0].size(); ii++){
      dfdu[0][ii] /= verts.size();
    }
  }
}

Eigen::MatrixXd FEM2DFun::dKdp(int eidx, int pidx)
{
  return Eigen::MatrixXd();
}

Eigen::VectorXd FEM2DFun::df()
{
  //gradient of f with respect to parameters.
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(param.size());
  //gradient of f with respect to material distribution.
  Eigen::VectorXd dfddist = Eigen::VectorXd::Zero(distribution.size());

  int nrows = (int)m_I.size() - 1;
  compute_dfdu();
  Eigen::MatrixXd K0 = (em->getStiffness(0)).cast<double>();
  //sensitivity analysis using the adjoint method.
  for (unsigned int ii = 0; ii < u.size(); ii++){
    std::vector<double> lambda(nrows, 0);
    //lambda = K^{-1} dfdu.
    //dfdx = -lambda * dK/dparam * u.
    std::vector<int> I = m_I;
    std::vector<int> J = m_J;
    //checkSparseIndex(I, J);
    sparseSolve(I.data(), J.data(), m_val.data(), nrows, &(lambda[0]), &(dfdu[ii][0]));
    for (unsigned int jj = 0; jj < em->e.size(); jj++){
      Element2D * ele = em->e[jj];
      int nV = ele->nV();
      Eigen::MatrixXd dKdp;
      Eigen::VectorXd Ue(nV * dim);
      Eigen::VectorXd lambda_e(nV * dim);
      //in this example dK/dParam = K * (1/(3-2*rho)^2).
      double numerator = (3 - 2 * distribution[jj]);
      numerator *= numerator;
      dKdp =  (3.0/numerator) * K0;
      for (int kk = 0; kk < ele->nV(); kk++){
        int vidx = em->e[jj]->at(kk);
        for (int ll = 0; ll < dim; ll++){
          Ue[kk*dim + ll] = u[ii][vidx*dim + ll];
          lambda_e[kk*dim + ll] = lambda[vidx*dim + ll];
        }
      }
      dfddist[jj] += -lambda_e.dot(dKdp * Ue);
    }
  }

  for (unsigned int ii = 0; ii < dfddist.size(); ii++){
    dfddist[ii] += (mw / distribution.size()) * (density - m0);
  }

  Eigen::VectorXd coord = Eigen::VectorXd::Zero(dim);
  for (int ii = 0; ii < m_nx; ii++){
    for (int jj = 0; jj < m_ny; jj++){
      int eidx = grid[ii][jj];
      coord[0] = (ii + 0.5) / m_nx;
      coord[1] = (jj + 0.5) / m_ny;
      Eigen::SparseVector<double> ddistdp = field->df(coord);
      for (Eigen::SparseVector<double>::InnerIterator it(ddistdp); it; ++it){
        grad[it.index()] += dfddist[eidx] * it.value();
      }
    }
  }
  return grad;
}

void FEM2DFun::log(std::ostream & out)
{
  out << dx << " " << dy << " " << density << "\n";
}

FEM2DFun::FEM2DFun() :em(0), dim(2),
m_Init(false),
m_periodic(true),
m_fixRigid(true),
dx0(1e-2), dy0(5e-3),
dxw(5), dyw(1),
forceMagnitude(100),
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
  MatrixXS K0 = em->getStiffness(0);
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    Element2D * ele = em->e[ii];
    int nV = ele->nV();
    //K = K0 * rho/(3-2*rho).
    //change for better bound if there is one
    MatrixXS K = (cfgScalar) (param[ii] / (3 - 2 * param[ii])) * K0;

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
            if (2 * vj + dim1 == 2 * vk + dim2){
              val *= 1 + 1e-5;
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
  int nx = (int)grid.size();
  int ny = (int)grid[0].size();
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
  int nx = (int)grid.size();
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
  int ny = (int)grid[0].size();
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
  int nx = (int)grid.size();
  int ny = (int)grid[0].size();
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
  Eigen::Vector3d fv = ff / (double)rightv.size();
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

void stretchY(ElementMesh2D * em, const Eigen::Vector3d & ff, const Grid2D& grid, std::vector<double> & externalForce)
{
  int dim = 2;
  std::vector<int> topv, botv;
  topv = topVerts(em, grid);
  Eigen::Vector3d fv = ff / (double)topv.size();
  for (unsigned int ii = 0; ii<topv.size(); ii++){
    int vidx = topv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx + jj] += fv[jj];
    }
  }
  botv = botVerts(em, grid);
  for (unsigned int ii = 0; ii<botv.size(); ii++){
    int vidx = botv[ii];
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

double measureShearX(ElementMesh2D * em, const std::vector<double> & u, const Grid2D & grid)
{
  int dim = 2;
  std::vector<int> tv, bv;
  tv = topVerts(em, grid);
  bv = botVerts(em, grid);
  double shear = 0;
  for (unsigned int ii = 0; ii<tv.size(); ii++){
    shear += u[dim * tv[ii]] - u[dim * bv[ii]];
  }
  shear /= bv.size();
  return shear;
}
