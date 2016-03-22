#include "ArrayUtil.hpp"
#include "FEM2DFun.hpp"
#include "Element2D.h"
#include "ElementMesh2D.h"
#include "RealField.hpp"
#include "pardiso_sym.hpp"
//#include "Timer.hpp"

typedef Eigen::Triplet<cfgScalar> TripletS;

///@brief a copy of stiffness matrix assembly function in ElementMesh2D class.
///Difference is that this function scales K for each element using param.
void getStiffnessSparse(ElementMesh2D * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic);

std::vector<int> topVerts(ElementMesh2D * em, const std::vector<int> & s);
std::vector<int> botVerts(ElementMesh2D * em, const std::vector<int> & s);
std::vector<int> leftVerts(ElementMesh2D * em, const std::vector<int> & s);
std::vector<int> rightVerts(ElementMesh2D * em, const std::vector<int> & s);
std::vector<int> cornerVerts(ElementMesh2D * em, const std::vector<int> & gridSize);
static const int sw2d[4][2] =
{ { -1, -1 },
{  -1, 1 },
{  1, -1 },
{  1, 1 }};
Eigen::MatrixXd BMatrix(const Vector2d & xx, const Eigen::Vector2d & size);
Eigen::VectorXd quadStrain(const Eigen::VectorXd & u, const Eigen::VectorXd & X,
  const Eigen::Vector2d & xi);
void copyVert2(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<Vector2S> & X);

  void copyVert2(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<double> & u);

void FEM2DFun::initArrays()
{
  int nForce = (int)externalForce.size();
  int nrow = (int)externalForce[0].size();
  u.resize(nForce);
  dfdu.resize(u.size());
  for (unsigned int ii = 0; ii < u.size(); ii++){
    u[ii].resize(nrow);
    dfdu[ii].resize(nrow);
  }
  G = Eigen::MatrixXd(3, nForce);
}

void FEM2DFun::init(const Eigen::VectorXd & x0)
{
  int nForce = 3;
  bool triangular = true;
  m_I.clear();
  m_J.clear();
  em->stiffnessPattern(m_I, m_J, triangular, m_fixRigid, m_fixRigid, m_periodic);
  //pardiso uses 1-based
  for (unsigned int ii = 0; ii < m_I.size(); ii++){
    m_I[ii] ++;
  }
  for (unsigned int ii = 0; ii < m_J.size(); ii++){
    m_J[ii] ++;
  }

  pardisoState = new PardisoState();
  pardisoInit(pardisoState);
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), (int)m_I.size() - 1, pardisoState);
  param = x0;
  distribution = Eigen::VectorXd::Zero(em->e.size());
  int nrow = (int)m_I.size() - 1;
  externalForce.resize(nForce);
  for (int ii = 0; ii < (int)externalForce.size(); ii++){
    externalForce[ii].resize(nrow, 0);
  }
  Eigen::Vector3d forceDir = forceMagnitude * Eigen::Vector3d(1, 0, 0);
  stretchX(em, forceDir, gridSize, externalForce[0]);
  forceDir = forceMagnitude * Eigen::Vector3d(0, 1, 0);
  stretchY(em, forceDir, gridSize, externalForce[1]);
  shearXY(em, forceMagnitude, gridSize, externalForce[2]);
  G0 = Eigen::MatrixXd::Identity(3, 3);
  wG = Eigen::VectorXd::Ones(3);
  initArrays();
  K0 = em->getStiffness(0);
  m_Init = true;
}

void FEM2DFun::setParam(const Eigen::VectorXd & x0)
{
  bool triangle = true;
  Vector3S color0(1, 1, 1);
  param = x0;

  //read material distribution from field
  for (int ii = 0; ii < field->param.size(); ii++){
    field->setParam(ii, x0[ii]);
  }
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      Eigen::VectorXd coord = Eigen::VectorXd::Zero(dim);
      int eIdx = gridToLinearIdx(ii, jj, gridSize);
      coord[0] = (ii + 0.5) / gridSize[0];
      coord[1] = (jj + 0.5) / gridSize[1];
      distribution[eIdx] = field->f(coord);
      em->e[eIdx]->color = distribution[eIdx] * color0;
    }
  }

  //solve linear statics problems
  //Timer timer;
  //timer.startWall();
  getStiffnessSparse();
  //timer.endWall();
  //std::cout << "assemble time " << timer.getSecondsWall() << "\n";
  int nrows = (int)m_I.size() - 1;

  pardisoNumericalFactorize(m_I.data(), m_J.data(), m_val.data(), nrows, pardisoState);
  //checkSparseIndex(m_I, m_J);
  for (unsigned int ii = 0; ii < externalForce.size(); ii++){
    std::fill(u[ii].begin(), u[ii].end(), 0);
    //timer.startWall();
    pardisoBackSubstitute(m_I.data(), m_J.data(), m_val.data(), nrows, u[ii].data(), externalForce[ii].data(), pardisoState);
    //timer.endWall();
    //std::cout<<"Lin solve time " << timer.getSecondsWall() << "\n";
  }
  
  //density objective
  density = 0;
  for (unsigned int ii = 0; ii < distribution.size(); ii++){
    density += distribution[ii];
  }
  density /= distribution.size();
  //std::cout << "dx dy: " << dx << " " << dy << "\n";
  std::vector<int> vidx = cornerVerts(em, gridSize);
  Eigen::VectorXd X;
  copyVert2(X, vidx, em->X);
  for (int ii = 0; ii < G.cols(); ii++){
    Eigen::Vector2d xi(0, 0);
    Eigen::VectorXd x;
    copyVert2(x, vidx, u[ii]);
    Eigen::VectorXd strain = quadStrain(x, X, xi);
    G.col(ii) = strain;
  }
  std::cout << G(0, 0) << " " << G(1, 0) << " " << G(2, 0) << "\n";
  std::cout << G(0, 1) << " " << G(1, 1) << " " << G(2, 1) << "\n";
  std::cout << G(0, 2) << " " << G(1, 2) << " " << G(2, 2) << "\n";
}

void copyVert2(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<Vector2S> & X)
{
  int dim = 2;
  x = Eigen::VectorXd(dim * vidx.size());
  for (unsigned int ii = 0; ii < vidx.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      x[ii*dim + jj] = X[vidx[ii]][jj];
    }
  }
}

void copyVert2(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<double> & u)
{
  int dim = 2;
  x = Eigen::VectorXd(dim * vidx.size());
  for (unsigned int ii = 0; ii < vidx.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      x[ii*dim + jj] = u[dim*vidx[ii] + jj];
    }
  }
}

double FEM2DFun::f()
{
  double val = 0;
  Eigen::MatrixXd diff = G - G0;
  diff = diff.cwiseProduct(diff);
  Eigen::VectorXd p = diff * wG;
  val = 0.5 * p.sum();
  val += 0.5 * mw * (density - m0) * (density - m0);
  //std::cout << val << "\n";
  return val;
}

void FEM2DFun::compute_dfdu()
{
  int nForce = G.cols();
  //only 8 corner vertices affect G
  std::vector<int> vidx = cornerVerts(em, gridSize);
  int nV = 4;
  int dim = 2;
  Eigen::VectorXd X(nV*dim);
  copyVert2(X, vidx, em->X);
  Eigen::Vector2d esize = X.segment<2>(2 * (nV - 1)) - X.segment<2>(0);
  Eigen::Vector2d xi(0, 0);
  Eigen::MatrixXd B = BMatrix(xi, esize);
  Eigen::MatrixXd diff = (G - G0) * wG.asDiagonal();
  Eigen::MatrixXd grad = B.transpose() * diff;
  //loop each displacement
  for (int ii = 0; ii < nForce; ii++){
    //each column corresponds to one set of u.
    std::fill(dfdu[ii].begin(), dfdu[ii].end(), 0);
    //loop 8 corners
    for (int jj = 0; jj < nV; jj++){
      for (int kk = 0; kk < dim; kk++){
        //each corner affects 6 entries in G.
        dfdu[ii][vidx[jj] * dim + kk] += grad(jj*dim + kk, ii);
      }
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
  //sensitivity analysis using the adjoint method.
  for (unsigned int ii = 0; ii < u.size(); ii++){
    std::vector<double> lambda(nrows, 0);
    //lambda = K^{-1} dfdu.
    //dfdx = -lambda * dK/dparam * u.
    pardisoBackSubstitute(m_I.data(), m_J.data(), m_val.data(), nrows, lambda.data(), dfdu[ii].data(), pardisoState);
    for (unsigned int jj = 0; jj < em->e.size(); jj++){
      Element2D * ele = em->e[jj];
      int nV = ele->nV();
      Eigen::MatrixXf dKdp;
      Eigen::VectorXf Ue(nV * dim);
      Eigen::VectorXf lambda_e(nV * dim);
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
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      int eidx = gridToLinearIdx(ii,jj,gridSize);
      coord[0] = (ii + 0.5) / gridSize[0];
      coord[1] = (jj + 0.5) / gridSize[1];
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
  out << G(0, 0) << " " << G(1, 0) << "\n";
}

FEM2DFun::FEM2DFun() :em(0), dim(2),
m_Init(false),
m_periodic(true),
m_fixRigid(true),
constrained(false),
forceMagnitude(1),
gridSize(2,0),
field(0)
{
}

FEM2DFun::~FEM2DFun(){}

MatrixXS FEM2DFun::getKe(int ei)
{
  return (cfgScalar)(distribution[ei]) * K0;
  //(cfgScalar)(distribution[ei] / (3 - 2 * distribution[ei])) * K0;
}

void FEM2DFun::getStiffnessSparse()
{
  bool trig = true;
  int N = 2 * (int)em->x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N, N);
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    Element2D * ele = em->e[ii];
    int nV = ele->nV();
    //K = K0 * rho/(3-2*rho).
    //change for better bound if there is one
    MatrixXS K = getKe(ii);

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
  if (m_fixRigid)
  {
    em->fixTranslation(Ksparse, trig, em);
    em->fixRotation(Ksparse, trig, em);
  }
  if (m_periodic)
  {
    em->enforcePeriodicity(Ksparse, trig, em);
  }

  m_val.resize(Ksparse.nonZeros());
  int idx = 0;
  for (int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
      m_val[idx] = it.value();
      idx++;
    }
  }
}

std::vector<int> topVerts(ElementMesh2D * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int topV[2] = { 1, 3 };
  for (int ii = 0; ii<nx; ii++){
    int ei = gridToLinearIdx(ii, ny - 1, s);
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(topV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

std::vector<int> botVerts(ElementMesh2D * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int botV[2] = { 0, 2 };
  for (int ii = 0; ii<nx; ii++){
    int ei = gridToLinearIdx(ii,0,s);
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(botV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

std::vector<int> leftVerts(ElementMesh2D * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int ny = s[1];
  int leftV[2] = { 0, 1 };
  for (int ii = 0; ii<ny; ii++){
    int ei = gridToLinearIdx(0,ii,s);
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(leftV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

std::vector<int> rightVerts(ElementMesh2D * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int rightV[2] = { 2, 3 };
  for (int ii = 0; ii<ny; ii++){
    int ei = gridToLinearIdx(nx - 1,ii,s);
    for (int jj = 0; jj<2; jj++){
      int vi = em->e[ei]->at(rightV[jj]);
      v.push_back(vi);
    }
  }
  return v;
}

void stretchX(ElementMesh2D * em, const Eigen::Vector3d & ff, const std::vector<int> & s, 
  std::vector<double> & externalForce)
{
  int dim = 2;
  std::vector<int> leftv, rightv;
  rightv = rightVerts(em, s);
  Eigen::Vector3d fv = ff / (double)rightv.size();
  for (unsigned int ii = 0; ii<rightv.size(); ii++){
    int vidx = rightv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx+jj] += fv[jj];
    }
  }
  leftv = leftVerts(em, s);
  for (unsigned int ii = 0; ii<leftv.size(); ii++){
    int vidx = leftv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx + jj] -= fv[jj];
    }
  }
}

void stretchY(ElementMesh2D * em, const Eigen::Vector3d & ff, const std::vector<int>& s, std::vector<double> & externalForce)
{
  int dim = 2;
  std::vector<int> topv, botv;
  topv = topVerts(em, s);
  Eigen::Vector3d fv = ff / (double)topv.size();
  for (unsigned int ii = 0; ii<topv.size(); ii++){
    int vidx = topv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx + jj] += fv[jj];
    }
  }
  botv = botVerts(em, s);
  for (unsigned int ii = 0; ii<botv.size(); ii++){
    int vidx = botv[ii];
    for (int jj = 0; jj < dim; jj++){
      externalForce[dim * vidx + jj] -= fv[jj];
    }
  }
}

void shearXY(ElementMesh2D * em, double ff, const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(1, 0, 0);
  //apply horizontal force on top and bottom faces
  stretchY(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 1, 0);
  stretchX(em, force, s, fe);
}

std::vector<int> cornerVerts(ElementMesh2D * em, const std::vector<int> & gridSize)
{
  std::vector<int> corner(4, 0);
  int x = gridSize[0] - 1;
  int y = gridSize[1] - 1;
  corner[0] = em->e[gridToLinearIdx( 0, 0, gridSize)]->at(0);
  corner[1] = em->e[gridToLinearIdx( 0, y, gridSize)]->at(1);
  corner[2] = em->e[gridToLinearIdx( x, 0, gridSize)]->at(2);
  corner[3] = em->e[gridToLinearIdx( x, y, gridSize)]->at(3);
  return corner;
}


//assuming element size 1.
Vector2d shapeFunGrad(int ii, const Vector2d & xx)
{
  Vector2d grad;
  grad[0] = sw2d[ii][0] * (1 + sw2d[ii][1] * xx[1]);
  grad[1] = sw2d[ii][1] * (1 + sw2d[ii][0] * xx[0]);
  return 0.5*grad;
}

Eigen::MatrixXd BMatrix(const Vector2d & xx, const Eigen::Vector2d & size)
{
  int nV = 4;
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * nV);
  for (int ii = 0; ii < nV; ii++){
    int col = 2 * ii;
    Vector2d dN = shapeFunGrad(ii, xx).cwiseQuotient(size);
    B(0, col) = dN[0];
    B(1, col + 1) = dN[1];
 
    B(2, col) = dN[1];
    B(2, col + 1) = dN[0];
  }
  return B;
}

//returns 6 element vector measured at xi.
Eigen::VectorXd quadStrain(const Eigen::VectorXd & u, const Eigen::VectorXd & X,
  const Eigen::Vector2d & xi)
{
  int dim = 2;
  //size of quad element
  int nV = X.rows() / 2;
  Eigen::Vector2d esize = X.segment<2>(2 * (nV - 1)) - X.segment<2>(0);
  Eigen::MatrixXd B = BMatrix(xi, esize);
  Eigen::VectorXd strain = B*u;
  return strain;
}
