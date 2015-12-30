#include "cfgDefs.h"
#include "ArrayUtil.hpp"
#include "FEM3DFun.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "RealField.hpp"
#include "pardiso_sym.hpp"
#include "Timer.hpp"

typedef Eigen::Triplet<cfgScalar> TripletS;

static const int sw[8][3] =
{ { -1, -1, -1 },
{ -1, -1, 1 },
{ -1, 1, -1 },
{ -1, 1, 1 },
{ 1, -1, -1 },
{ 1, -1, 1 },
{ 1, 1, -1 },
{ 1, 1, 1 }
};

///@brief a copy of stiffness matrix assembly function in ElementMesh2D class.
///Difference is that this function scales K for each element using param.
void getStiffnessSparse(ElementMesh * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic);

std::vector<int> topVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> botVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> leftVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> rightVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> frontVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> backVerts(ElementMesh * em, const std::vector<int> & gridSize);

///@brief add f to each 3 subvector of a.
void addVector3d(std::vector<double> & a, const Eigen::Vector3d & f,
  const std::vector<int> & idx);

///@brief eight corners of a grid for coarsening.
std::vector<int> cornerVerts(ElementMesh * em, const std::vector<int> & gridSize);

///@brief assuming element size 1.
Vector3d shapeFunGrad(int ii, const Vector3d & xx);

///@param size. size of a cube element.
Eigen::MatrixXd BMatrix(const Vector3d & xx, const Eigen::Vector3d & size);

int gridToLinearIdx(int ix, int iy, int iz, const std::vector<int> & gridSize)
{
  return ix * gridSize[1] * gridSize[2] + iy * gridSize[2] + iz;
}

void FEM3DFun::initArrays()
{
  int nForce = (int)externalForce.size();
  int nrow = (int)externalForce[0].size();
  u.resize(nForce);
  dfdu.resize(u.size());
  for (unsigned int ii = 0; ii < u.size(); ii++){
    u[ii].resize(nrow);
    dfdu[ii].resize(nrow);
  }
  G = Eigen::MatrixXd(6, nForce);
}

void FEM3DFun::init(const Eigen::VectorXd & x0)
{  
  //6 harmonic displacements.
  int nForce = 6;
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
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), m_I.size()-1, pardisoState);
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
  forceDir = forceMagnitude * Eigen::Vector3d(0, 0, 1);
  stretchZ(em, forceDir, gridSize, externalForce[2]);
  shearYZ(em, forceMagnitude, gridSize, externalForce[3]);
  shearXZ(em, forceMagnitude, gridSize, externalForce[4]);
  shearXY(em, forceMagnitude, gridSize, externalForce[5]);

  initArrays();
  m_Init = true;
}

void FEM3DFun::setParam(const Eigen::VectorXd & x0)
{
  bool triangle = true;
  bool constrained = false;
  Vector3S color0(0.6, 0.6, 1.0);
  param = x0;

  //read material distribution from field
  for (int ii = 0; ii < field->param.size(); ii++){
    field->setParam(ii, x0[ii]);
  }
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      for (int kk = 0; kk < gridSize[2]; kk++){
        Eigen::VectorXd coord = Eigen::VectorXd::Zero(dim);
        int eIdx = gridToLinearIdx(ii, jj, kk, gridSize);
        coord[0] = (ii + 0.5) / gridSize[0];
        coord[1] = (jj + 0.5) / gridSize[1];
        coord[2] = (kk + 0.5) / gridSize[2];
        distribution[eIdx] = field->f(coord);
        em->e[eIdx]->color = distribution[eIdx] * color0;
      }
    }
  }

  //solve linear statics problems
  Timer timer;
  std::vector<cfgScalar> val;
  //timer.startWall();
  getStiffnessSparse(em, distribution, val, triangle, constrained, m_fixRigid, m_periodic);
  //timer.endWall();
  //std::cout << "assemble time " << timer.getSecondsWall() << "\n";
  m_val = std::vector<double>(val.begin(), val.end());
  int nrows = (int)m_I.size() - 1;

  timer.startWall();
  pardisoNumericalFactorize(m_I.data(), m_J.data(), m_val.data(), nrows, pardisoState);
  timer.endWall();
  std::cout << "num fact time " << timer.getSecondsWall() << "\n";

  for (unsigned int ii = 0; ii < externalForce.size(); ii++){
    std::fill(u[ii].begin(), u[ii].end(), 0);
    timer.startWall();
    //checkSparseIndex(I, J);
    pardisoBackSubstitute(m_I.data(), m_J.data(), m_val.data(), nrows, u[ii].data(), externalForce[ii].data(), pardisoState);
    timer.endWall();
    std::cout<<"Lin subst time " << timer.getSecondsWall() << "\n";
  }

  //density objective
  density = 0;
  for (unsigned int ii = 0; ii < distribution.size(); ii++){
    density += distribution[ii];
  }
  density /= distribution.size();

}

double FEM3DFun::f()
{
  double val = 0;
  return val ;
}

void FEM3DFun::compute_dfdu()
{
  
}

Eigen::MatrixXd FEM3DFun::dKdp(int eidx, int pidx)
{
  return Eigen::MatrixXd();
}

Eigen::VectorXd FEM3DFun::df()
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
    pardisoBackSubstitute(m_I.data(), m_J.data(), m_val.data(), nrows, lambda.data(), dfdu[ii].data(), pardisoState);
    for (unsigned int jj = 0; jj < em->e.size(); jj++){
      Element * ele = em->e[jj];
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
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      for (int kk = 0; kk < gridSize[2]; kk++){
        int eidx = gridToLinearIdx(ii, jj, kk, gridSize);
        coord[0] = (ii + 0.5) / gridSize[0];
        coord[1] = (jj + 0.5) / gridSize[1];
        coord[2] = (kk + 0.5) / gridSize[2];
        Eigen::SparseVector<double> ddistdp = field->df(coord);
        for (Eigen::SparseVector<double>::InnerIterator it(ddistdp); it; ++it){
          grad[it.index()] += dfddist[eidx] * it.value();
        }
      }
    }
  }
  return grad;
}

void FEM3DFun::log(std::ostream & out)
{
  
}

FEM3DFun::FEM3DFun() :em(0), dim(3),
m_Init(false),
m_periodic(true),
m_fixRigid(true),
forceMagnitude(100),
gridSize(3,0),
field(0)
{
}

FEM3DFun::~FEM3DFun(){}

void getStiffnessSparse(ElementMesh * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic)
{
  int dim = 3;
  int N = dim * (int)em->x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N, N);
  MatrixXS K0 = em->getStiffness(0);
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    Element * ele = em->e[ii];
    int nV = ele->nV();
    //K = K0 * rho/(3-2*rho).
    //change for better bound if there is one
    MatrixXS K = (cfgScalar) (param[ii] / (3 - 2 * param[ii])) * K0;

    for (int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for (int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for (int dim1 = 0; dim1<dim; dim1++){
          for (int dim2 = 0; dim2<dim; dim2++){
            if (trig && (dim * vk + dim2 > dim * vj + dim1)) {
              continue;
            }
            cfgScalar val = K(dim * jj + dim1, dim * kk + dim2);
            if (constrained){
              if ( (em->fixed[vk] || em->fixed[vj]) 
                && (vj != vk || dim1 != dim2)){
                val = 0;
              }
            }
            if (dim * vj + dim1 == dim * vk + dim2){
              val *= 1 + 1e-5;
            }
            TripletS triple(dim * vj + dim1, dim * vk + dim2, val);
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

  val.resize(Ksparse.nonZeros());
  int idx = 0;
  for (int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
      val[idx] = it.value();
      idx++;
    }
  }
}

std::vector<int> topVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int nz = s[2];
  int topV[4] = { 2, 3, 6, 7 };
  for (int ii = 0; ii<nx; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(ii, ny - 1, kk, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(topV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> botVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int nz = s[2];
  int botV[4] = { 0, 1, 4, 5};
  for (int ii = 0; ii<nx; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(ii,0,kk,s);
      for (int jj = 0; jj<4; jj++){
        int vi = em->e[ei]->at(botV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> leftVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int ny = s[1];
  int nz = s[2];
  int leftV[4] = { 0, 1 ,2 ,3};
  for (int ii = 0; ii<ny; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(0,ii,kk, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(leftV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> rightVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int nz = s[2];
  int rightV[4] = { 4, 5, 6, 7 };
  for (int ii = 0; ii<ny; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(nx-1,ii,kk,s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(rightV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> frontVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int nz = s[2];
  int rightV[4] = { 1, 3, 5, 7 };
  for (int ii = 0; ii<nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      int ei = gridToLinearIdx(ii, jj, nz-1, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(rightV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> backVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int rightV[4] = { 0, 2, 4, 6 };
  for (int ii = 0; ii<nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      int ei = gridToLinearIdx(ii, jj, 0, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(rightV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

void addVector3d(std::vector<double> & a, const Eigen::Vector3d & f,
  const std::vector<int> & idx)
{
  int dim = f.size(); 
  for (unsigned int ii = 0; ii < idx.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      a[dim * idx[ii] + jj] += f[jj];
    }
  }
}

void stretchX(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int> & s, std::vector<double> & externalForce)
{
  std::vector<int> leftv, rightv;
  rightv = rightVerts(em, s);
  Eigen::Vector3d fv = ff / (double)rightv.size();
  addVector3d(externalForce, fv, rightv);
  leftv = leftVerts(em, s);
  addVector3d(externalForce, -fv, leftv);
}

void stretchY(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int> & s, std::vector<double> & externalForce)
{
  std::vector<int> topv, botv;
  topv = topVerts(em, s);
  Eigen::Vector3d fv = ff / (double)topv.size();
  addVector3d(externalForce, fv, topv);
  botv = botVerts(em, s);
  addVector3d(externalForce, -fv, botv);
}

void stretchZ(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int>& s, std::vector<double> & externalForce)
{
  std::vector<int> frontv, backv;
  frontv = frontVerts(em, s);
  Eigen::Vector3d fv = ff / (double)frontv.size();
  addVector3d(externalForce, fv, frontv);
  backv = backVerts(em, s);
  addVector3d(externalForce, -fv, backv);
}

void shearXY(ElementMesh * em, double ff,
    const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(1, 0, 0);
  //apply horizontal force on top and bottom faces
  stretchY(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 1, 0);
  stretchX(em, force, s, fe);
}

void shearYZ(ElementMesh * em, double ff,
  const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(0, 1, 0);
  //apply horizontal force on top and bottom faces
  stretchZ(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 0, 1);
  stretchY(em, force, s, fe);
}

void shearXZ(ElementMesh * em, double ff, 
  const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(1, 0, 0);
  //apply horizontal force on top and bottom faces
  stretchZ(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 0, 1);
  stretchX(em, force, s, fe);
}

std::vector<int> cornerVerts(ElementMesh * em, const std::vector<int> & gridSize)
{
  std::vector<int> corner(8, 0);
  int x = gridSize[0] - 1;
  int y = gridSize[1] - 1;
  int z = gridSize[2] - 1;
  corner[0] = em->e[gridToLinearIdx(0, 0, 0, gridSize)]->at(0);
  corner[1] = em->e[gridToLinearIdx(0, 0, z, gridSize)]->at(1);
  corner[2] = em->e[gridToLinearIdx(0, y, 0, gridSize)]->at(2);
  corner[3] = em->e[gridToLinearIdx(0, y, z, gridSize)]->at(3);
  corner[4] = em->e[gridToLinearIdx(x, 0, 0, gridSize)]->at(4);
  corner[5] = em->e[gridToLinearIdx(x, 0, z, gridSize)]->at(5);
  corner[6] = em->e[gridToLinearIdx(x, y, 0, gridSize)]->at(6);
  corner[7] = em->e[gridToLinearIdx(x, y, z, gridSize)]->at(7);
  return corner;
}

//assuming element size 1.
Vector3d shapeFunGrad(int ii, const Vector3d & xx)
{
  Vector3d grad;
  grad[0] = sw[ii][0] * (1 + sw[ii][1] * xx[1]) * (1 + sw[ii][2] * xx[2]);
  grad[1] = sw[ii][1] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][2] * xx[2]);
  grad[2] = sw[ii][2] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][1] * xx[1]);
  return 0.25*grad;
}

Eigen::MatrixXd BMatrix(const Vector3d & xx, const Eigen::Vector3d & size)
{
  int nV = 8;
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 3 * nV);
  for (int ii = 0; ii < nV; ii++){
    int col = 3 * ii;
    Vector3d dN = shapeFunGrad(ii, xx).cwiseQuotient(size);
    B(0, col) = dN[0];
    B(1, col + 1) = dN[1];
    B(2, col + 2) = dN[2];

    B(3, col) = dN[1];
    B(3, col + 1) = dN[0];

    B(4, col + 1) = dN[2];
    B(4, col + 2) = dN[1];

    B(5, col) = dN[2];
    B(5, col + 2) = dN[0];
  }
  return B;
}

//returns 6 element vector measured at xi.
Eigen::VectorXd hexStrain(const Eigen::VectorXd & x, const Eigen::VectorXd & X,
  const Eigen::Vector3d & xi)
{
  int dim = 3;
  //size of cube element
  int nV = X.cols()/3;
  Eigen::Vector3d esize = X.segment<3>(0) - X.segment<3>(3*(nV-1));
  Eigen::MatrixXd B=BMatrix(xi, esize);
  Eigen::VectorXd strain = B*(x-X);
  return strain;
}
