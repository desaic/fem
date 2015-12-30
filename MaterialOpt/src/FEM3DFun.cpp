#include "cfgDefs.h"
#include "ArrayUtil.hpp"
#include "FEM3DFun.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "RealField.hpp"
#include "SparseLin.hpp"
//#include "Timer.hpp"

typedef Eigen::Triplet<cfgScalar> TripletS;

///@brief a copy of stiffness matrix assembly function in ElementMesh2D class.
///Difference is that this function scales K for each element using param.
void getStiffnessSparse(ElementMesh * em, const Eigen::VectorXd & param,
  std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid, bool iPeriodic);

std::vector<int> topVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> botVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> leftVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> rightVerts(ElementMesh * em, const std::vector<int> & gridSize);

int gridToLinearIdx(int ix, int iy, int iz, const std::vector<int> & gridSize)
{
  return ix * gridSize[1] * gridSize[2] + iy * gridSize[2] + iz;
}

void FEM3DFun::initArrays()
{
  int nrow = externalForce[0].size();
  u.resize(externalForce.size());
  dfdu.resize(u.size());
  for (unsigned int ii = 0; ii < u.size(); ii++){
    u[ii].resize(nrow);
    dfdu[ii].resize(nrow);
  }
}

void FEM3DFun::init(const Eigen::VectorXd & x0)
{  
  bool triangular = true;
  m_I.clear();
  m_J.clear();
  em->stiffnessPattern(m_I, m_J, triangular, m_fixRigid, m_fixRigid, m_periodic);
  sparseInit();
  param = x0;
  distribution = Eigen::VectorXd::Zero(em->e.size());
  int nrow = (int)m_I.size() - 1;
  externalForce.resize(1);
  externalForce[0].resize(nrow, 0);
  Eigen::Vector3d forceDir = forceMagnitude * Eigen::Vector3d(1, 0, 0);
  stretchX(em, forceDir, gridSize, externalForce[0]);
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
  //Timer timer;
  std::vector<cfgScalar> val;
  //timer.startWall();
  getStiffnessSparse(em, distribution, val, triangle, constrained, m_fixRigid, m_periodic);
  //timer.endWall();
  //std::cout << "assemble time " << timer.getSecondsWall() << "\n";
  m_val = std::vector<double>(val.begin(), val.end());
  int nrows = (int)m_I.size() - 1;
  
  for (unsigned int ii = 0; ii < externalForce.size(); ii++){
    std::fill(u[ii].begin(), u[ii].end(), 0);
    std::vector<int> I = m_I;
    std::vector<int> J = m_J;
    //timer.startWall();
    //checkSparseIndex(I, J);
    sparseSolve(I.data(), J.data(), m_val.data(), nrows, &(u[ii][0]), externalForce[ii].data());
    //timer.endWall();
    //std::cout<<"Lin solve time " << timer.getSecondsWall() << "\n";
  }

  dx = measureStretchX(em, u[0], gridSize);
  dy = measureStretchY(em, u[0], gridSize);

  //show rendering
  for (unsigned int ii = 0; ii < em->x.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      em->x[ii][jj] = em->X[ii][jj] + u[0][ii*dim + jj];
      em->fe[ii][jj] = externalForce[0][ii*dim + jj];
    }
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
  //Example: measure width and height under a stretching force.
  double val = 0.5 * dxw * (dx - dx0) * (dx - dx0) + 0.5 * dyw * (dy - dy0) * (dy - dy0);
  val += 0.5 * mw * (density - m0) * (density - m0);
  return val ;
}

void FEM3DFun::compute_dfdu()
{
  int nrows = (int)m_I.size() - 1;
  std::fill(dfdu[0].begin(), dfdu[0].end(), 0);
  for (unsigned int ii = 0; ii < dfdu.size(); ii++){
    std::vector<int> verts;
    verts = topVerts(em, gridSize);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii] + 1] += (dy - dy0)*dyw;
    }
    verts = botVerts(em, gridSize);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii] + 1] -= (dy - dy0)*dyw;
    }
    verts = rightVerts(em, gridSize);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii]] += (dx - dx0)*dxw;
    }
    verts = leftVerts(em, gridSize);
    for (unsigned int ii = 0; ii<verts.size(); ii++){
      dfdu[0][dim * verts[ii]] -= (dx - dx0)*dxw;
    }
    for (unsigned int ii = 0; ii<dfdu[0].size(); ii++){
      dfdu[0][ii] /= verts.size();
    }
  }
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
    sparseSolve(I.data(), J.data(), m_val.data(), nrows, &(lambda[0]), &(dfdu[ii][0]));
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
  out << dx << " " << dy << " " << density << "\n";
}

FEM3DFun::FEM3DFun() :em(0), dim(3),
m_Init(false),
m_periodic(true),
m_fixRigid(true),
dx0(1e-2), dy0(5e-3),
dxw(5), dyw(1),
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

  for (int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
      val.push_back(it.value());
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

void stretchX(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int> & s, std::vector<double> & externalForce)
{
  int dim = 3;
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

void stretchY(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int> & s, std::vector<double> & externalForce)
{
  int dim = 3;
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

double measureStretchX(ElementMesh * em, const std::vector<double> & u, const std::vector<int> & s)
{
  int dim = 3;
  std::vector<int> lv, rv;
  lv = leftVerts(em, s);
  rv = rightVerts(em, s);
  double stretch = 0;
  for (unsigned int ii = 0; ii<lv.size(); ii++){
    stretch += u[dim * rv[ii]] - u[dim * lv[ii]];
  }
  stretch /= lv.size();
  return stretch;
}

double measureStretchY(ElementMesh * em, const std::vector<double> & u, const std::vector<int> & s)
{
  int dim = 3;
  std::vector<int> tv, bv;
  tv = topVerts(em, s);
  bv = botVerts(em, s);
  double stretch = 0;
  for (unsigned int ii = 0; ii<tv.size(); ii++){
    stretch += u[dim * tv[ii]+1] - u[dim * bv[ii]+1];
  }
  stretch /= bv.size();
  return stretch;
}

double measureShearXY(ElementMesh * em, const std::vector<double> & u, const std::vector<int> & s)
{
  int dim = 3;
  std::vector<int> tv, bv;
  tv = topVerts(em, s);
  bv = botVerts(em, s);
  double shear = 0;
  for (unsigned int ii = 0; ii<tv.size(); ii++){
    shear += u[dim * tv[ii]] - u[dim * bv[ii]];
  }
  shear /= bv.size();
  return shear;
}
