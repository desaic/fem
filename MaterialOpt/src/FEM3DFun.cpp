#include "cfgDefs.h"
#include "ArrayUtil.hpp"
#include "FEM3DFun.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "RealField.hpp"
#include "pardiso_sym.hpp"
#include "Timer.hpp"

#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_builder.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/coarsening/plain_aggregates.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/pointwise_aggregates.hpp>
#include <amgcl/relaxation/multicolor_gauss_seidel.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/profiler.hpp>
namespace amgcl {
  profiler<> prof;
}

typedef Eigen::Triplet<cfgScalar> TripletS;
void amgLinSolve(std::vector<int> & I, std::vector<int> &J, std::vector<double> & val, std::vector<double> & b,
  std::vector<double> & x);
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
  //triangular = false;
  triangular = true;
  //6 harmonic displacements.
  int nForce = 6;
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
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), (int)m_I.size()-1, pardisoState);
  param = x0;
  distribution[0] = Eigen::VectorXd::Zero(em->e.size());
  distribution[1] = Eigen::VectorXd::Zero(em->e.size());
  int nrow = (int)m_I.size() - 1;
  externalForce.resize(nForce);
  for (int ii = 0; ii < (int)externalForce.size(); ii++){
    externalForce[ii].resize(nrow, 0);
  }
  
  //std::fill(externalForce[0].begin(), externalForce[0].end(), 0);
  //for (unsigned int ii = 0; ii < em->e.size(); ii++){
  //  for (int jj = 0; jj < em->e[ii]->nV(); jj++){
  //    int vidx = em->e[ii]->at(jj);
  //    for (int kk = 0; kk < dim; kk++){
  //      externalForce[0][vidx * dim + kk] += 1e-2 * forceMagnitude * twistY[jj][kk];
  //    }
  //  }
  //}
  Eigen::Vector3d forceDir = forceMagnitude * Eigen::Vector3d(1, 0, 0);
  stretchX(em, forceDir, gridSize, externalForce[0]);
  forceDir = forceMagnitude * Eigen::Vector3d(0, 1, 0);
  stretchY(em, forceDir, gridSize, externalForce[1]);
  forceDir = forceMagnitude * Eigen::Vector3d(0, 0, 1);
  stretchZ(em, forceDir, gridSize, externalForce[2]);
  shearXY(em, forceMagnitude, gridSize, externalForce[3]);
  shearYZ(em, forceMagnitude, gridSize, externalForce[4]);
  shearXZ(em, forceMagnitude, gridSize, externalForce[5]);
  G0 = Eigen::MatrixXd::Identity(6, 6);
  wG = Eigen::VectorXd::Ones(6);
  initArrays();
  K0 = em->getStiffness(0);
  m_Init = true;
}

void FEM3DFun::setParam(const Eigen::VectorXd & x0)
{
  Vector3S color0(0.8f, 0.8f, 1.0f);
  param = x0;

  //read material distribution from field
  for (int ii = 0; ii < field[0]->param.size(); ii++){
    field[0]->setParam(ii, x0[2*ii]);
    field[1]->setParam(ii, x0[2 * ii + 1]);
    //std::cout << x0[ii] << "\n";
  }
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      for (int kk = 0; kk < gridSize[2]; kk++){
        Eigen::VectorXd coord = Eigen::VectorXd::Zero(dim);
        int eIdx = gridToLinearIdx(ii, jj, kk, gridSize);
        coord[0] = (ii + 0.5) / gridSize[0];
        coord[1] = (jj + 0.5) / gridSize[1];
        coord[2] = (kk + 0.5) / gridSize[2];
        distribution[0][eIdx] = field[0]->f(coord);
        distribution[1][eIdx] = 1;// field[1]->f(coord);
        em->e[eIdx]->color = distribution[0][eIdx] * color0;
        em->e[eIdx]->color[2] *= distribution[1][eIdx];
        if (jj == 0){
          int b[4] = { 0, 1, 4, 5 };
          for (int ll = 0; ll < 4; ll++){
            em->fixed[em->e[eIdx]->at(b[ll])] = 1;
          }
        }
      }
    }
  }

  //solve linear statics problems
  Timer timer;
  timer.startWall();
  getStiffnessSparse();
  timer.endWall();
  std::cout << "assemble time " << timer.getSecondsWall() << "\n";
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
    //amgLinSolve(m_I, m_J, m_val, externalForce[ii], u[ii]);
    timer.endWall();
    //std::cout<<"Lin subst time " << timer.getSecondsWall() << "\n";
  }

  //density objective
  density = 0;
  for (unsigned int ii = 0; ii < distribution[0].size(); ii++){
    density += distribution[0][ii];
  }
  density /= distribution[0].size();

  std::vector<int> vidx = cornerVerts(em, gridSize);
  Eigen::VectorXd X;
  copyVert3(X, vidx, em->X);
  for (int ii = 0; ii < G.cols(); ii++){
    Eigen::Vector3d xi(0, 0, 0);
    Eigen::VectorXd x;
    copyVert3(x, vidx, u[ii]);
    Eigen::VectorXd strain = hexStrain(x, X, xi);
    G.col(ii) = strain;
  }

}

double FEM3DFun::f()
{
  double val = 0;
  Eigen::MatrixXd diff = G - G0;
  diff = diff.cwiseProduct(diff);
  Eigen::VectorXd p = diff * wG;
  val = 0.5 * p.sum();
  val += 0.5 * mw * (density - m0) * (density - m0);
  std::cout << G(0,0)<<" "<<G(1,0)<<" "<<G(2,0) << " "<<G(3,3)<<"\n";
  std::cout << val << "\n";
  return val ;
}

void FEM3DFun::compute_dfdu()
{
  int nForce = G.cols();
  //only 8 corner vertices affect G
  std::vector<int> vidx = cornerVerts(em, gridSize);
  int nV = 8;
  int dim = 3;
  Eigen::VectorXd X(nV*dim);
  copyVert3(X, vidx, em->X);
  Eigen::Vector3d esize = X.segment<3>(3 * (nV - 1)) - X.segment<3>(0);
  Eigen::Vector3d xi(0, 0, 0);
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

Eigen::MatrixXd FEM3DFun::dKdp(int eidx, int pidx)
{
  return Eigen::MatrixXd();
}

double dx_rho(double rho)
{
  //2 materioals
  return (1 - 1e-3) * 3 / ((3 - 2 * rho)*(3 - 2 * rho));
}

void dx_rho(double r0, double r1, double * drho)
{
  //3 materioals
  double p = 3;
  double E1 = 1;
  double E2 = 1;

  double x0 = std::pow(r0, p);
  double x01 = std::pow(r0, p-1);
  double x1 = 1;// std::pow(r1, p);
  double x11 = 0;// std::pow(r1, p - 1);
  double s;
  drho[0] = p * x01 * (x1*E1 + (1 - x1)*E2);
  s = p*x0 * (x11*E1 - x11*E2);
  drho[1] = 0;// s;

}

MatrixXS FEM3DFun::getKe(int ei)
{
  //return (cfgScalar)(distribution[ei] / (3 - 2 * distribution[ei])) * K0;
  //3 materials
  double p = 3;
  double E1 = 1;
  double E2 = 1;

  double x0 = 1e-3 + (1 - 1e-3)*std::pow(distribution[0][ei], p);
  double x1 = 1;// std::pow(distribution[1][ei], p);
  double s;
  s = x0;// *(x1*E1 + (1 - x1)*E2);
  return s*K0;
}

Eigen::VectorXd FEM3DFun::df()
{
  //gradient of f with respect to parameters.
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(param.size());
  //gradient of f with respect to material distribution.
  Eigen::VectorXd dfddist = Eigen::VectorXd::Zero(2*distribution[0].size());

  int nrows = (int)m_I.size() - 1;
  compute_dfdu();
  //sensitivity analysis using the adjoint method.
  for (unsigned int ii = 0; ii < u.size(); ii++){
    std::vector<double> lambda(nrows, 0);
    //lambda = K^{-1} dfdu.
    //dfdx = -lambda * dK/dparam * u.
    //checkSparseIndex(I, J);
    pardisoBackSubstitute(m_I.data(), m_J.data(), m_val.data(), nrows, lambda.data(), dfdu[ii].data(), pardisoState);
    for (unsigned int jj = 0; jj < em->e.size(); jj++){
      Element * ele = em->e[jj];
      int nV = ele->nV();
      Eigen::MatrixXf dKdp[2];
      Eigen::VectorXf Ue(nV * dim);
      Eigen::VectorXf lambda_e(nV * dim);
      double drho[2];
      dx_rho(distribution[0][ii], distribution[1][ii], drho);
      dKdp[0] =  drho[0] * K0;
      dKdp[1] = drho[1] * K0;
      
      for (int kk = 0; kk < ele->nV(); kk++){
        int vidx = em->e[jj]->at(kk);
        for (int ll = 0; ll < dim; ll++){
          Ue[kk*dim + ll] = u[ii][vidx*dim + ll];
          lambda_e[kk*dim + ll] = lambda[vidx*dim + ll];
        }
      }
      dfddist[2*jj] += -lambda_e.dot(dKdp[0] * Ue);
      dfddist[2 * jj+1] += -lambda_e.dot(dKdp[0] * Ue);
    }
  }

  for (unsigned int ii = 0; ii < distribution[0].size(); ii++){
    dfddist[2*ii] += (mw / distribution[0].size()) * (density - m0);
  }

  Eigen::VectorXd coord = Eigen::VectorXd::Zero(dim);
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      for (int kk = 0; kk < gridSize[2]; kk++){
        int eidx = gridToLinearIdx(ii, jj, kk, gridSize);
        coord[0] = (ii + 0.5) / gridSize[0];
        coord[1] = (jj + 0.5) / gridSize[1];
        coord[2] = (kk + 0.5) / gridSize[2];
        Eigen::SparseVector<double> ddistdp = field[0]->df(coord);
        for (Eigen::SparseVector<double>::InnerIterator it(ddistdp); it; ++it){
          grad[2* it.index()] += dfddist[2*eidx] * it.value();
          grad[2 * it.index() + 1] += dfddist[2 * eidx + 1] * it.value();
        }
      }
    }
  }
  return grad;
}

void FEM3DFun::log(std::ostream & out)
{
  out << G << "\n";
}

FEM3DFun::FEM3DFun() :em(0), dim(3),
m_Init(false),
m_periodic(true),
m_fixRigid(true),
constrained(false),
triangular(true),
forceMagnitude(0.5),
gridSize(3,0)
//field(0)
{
}

FEM3DFun::~FEM3DFun(){}


void FEM3DFun::getStiffnessSparse()
{
  bool trig = triangular;
  int dim = 3;
  int N = dim * (int)em->x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N, N);
  for (unsigned int ii = 0; ii<em->e.size(); ii++){
    Element * ele = em->e[ii];
    int nV = ele->nV();
    //K = K0 * rho/(3-2*rho).
    //change for better bound if there is one
    MatrixXS K = getKe(ii);

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

void amgLinSolve(std::vector<int> & I, std::vector<int> &J, std::vector<double> & val, std::vector<double> & b,
  std::vector<double> & x)
{
  using amgcl::prof;
  int n = I.size() - 1;
  
  prof.tic("build");
  typedef amgcl::make_solver<
    amgcl::amg<
    amgcl::backend::builtin<double>,
    amgcl::coarsening::aggregation,
    amgcl::relaxation::gauss_seidel
    >,
    amgcl::solver::cg<
    amgcl::backend::builtin<double>
    >
  > Solver;

  std::vector<int> ptr = I;
  std::vector<int> col = J;
  for (unsigned int ii = 0; ii < ptr.size(); ii++){
    ptr[ii] --;
  }
  for (unsigned int ii = 0; ii < col.size(); ii++){
    col[ii] --;
  }

  Solver::params prm;
  prm.precond.npre = 1;
  prm.precond.npost = 1;
  prm.precond.ncycle = 1;
  prm.solver.maxiter = 100;
  prm.precond.coarsening.aggr.block_size= 3;	////for pointwise aggregation
  prm.solver.tol = 1;
  prm.precond.coarse_enough = 2500;
  Solver solve(boost::tie(n, ptr, col, val), prm);
  prof.toc("build");

  std::cout << solve.precond() << std::endl;

  prof.tic("solve");
  size_t iters;
  double resid;
  boost::fill(x, 0);
  boost::tie(iters, resid) = solve(b, x);
  prof.toc("solve");

  std::cout << "Solver:" << std::endl
    << "  Iterations: " << iters << std::endl
    << "  Error:      " << resid << std::endl
    << std::endl;

  // Use the constructed solver as a preconditioner for another iterative
  // solver.
  //
  // Iterative methods use estimated residual for exit condition. For some
  // problems the value of estimated residual can get too far from true
  // residual due to round-off errors.
  //
  // Nesting iterative solvers in this way allows to shave last bits off the
  // error.
  //amgcl::solver::cg< amgcl::backend::builtin<double> > S(n);
  //boost::fill(x, 0);

  //prof.tic("nested solver");
  //boost::tie(iters, resid) = S(solve.system_matrix(), solve, b, x);
  //prof.toc("nested solver");

  //std::cout << "Nested solver:" << std::endl
  //  << "  Iterations: " << iters << std::endl
  //  << "  Error:      " << resid << std::endl
  //  << std::endl;

  std::cout << prof << std::endl;

}
