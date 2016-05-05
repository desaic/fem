#include "cfgDefs.h"
#include "ArrayUtil.hpp"
#include "FEM3DFun.hpp"
#include "ElementMeshUtil.hpp"
#include "RealField.hpp"
static const int d = 3;
Type_Define_VectorD_Outside_Class(d); Type_Define_VectorDi_Outside_Class(d); Type_Define_MatrixD_Outside_Class(d);
void FEM3DFun::initArrays()
{
  G = Eigen::MatrixXd(6, 6);
  G0 = G;
  wG = Eigen::MatrixXd::Zero(6, 6);
}

void FEM3DFun::init(const Eigen::VectorXd & x0)
{  
  m_Init = true;
}

double femCoef(double x){
  return x;// (x) / (3 - 2 * x);
}

double dfemCoef(double x){
  double denom = 3 - 2 * x;
  return 1;// 3 / (denom*denom);
}

///@param N grid size of u. 1 greater than u0.
///@param dx cell size. Assuming square cells.
void copyPeriodicU(int N, double dx, Eigen::VectorXd & u, const Eigen::VectorXd & u0,
  const Eigen::Matrix3d & strain)
{
  int dim = 3;
  for (int i = 0; i < N; i++){
    int i0 = i % (N - 1);
    for (int j = 0; j < N; j++){
      int j0 = j % (N - 1);
      for (int k = 0; k < N; k++){
        int k0 = k % (N - 1);
        Eigen::Vector3d X(i*dx, j*dx, k*dx);
        Eigen::Vector3d X0(i0*dx, j0*dx, k0*dx);
        Eigen::Vector3d y = strain * X;
        Eigen::Vector3d y0 = strain * X0;
        for (int d = 0; d < dim; d++){
          int idx = 3 * (i*N*N + j*N + k) + d;
          int idx0 = 3 * (i0 * (N-1)*(N-1) + j0 *(N-1) + k0) + d;
          u[idx] = u0[idx0];
          if ((i == N - 1) || (j == N - 1) || (k == N - 1)){
            u[idx] = u0[idx0] - y0[d] + y[d];
          }
        }
      }
    }
  }
}

static const int oneV[8][3] = { { 0, 0, 0 }, { 0, 0, 1 },
{ 0, 1, 0 }, { 0, 1, 1 }, { 1, 0, 0 }, { 1, 0, 1 },
{ 1, 1, 0 }, { 1, 1, 1 } };

std::vector<int> cell_verts(int cellIdx, int N)
{
  int nVert = 8;
  std::vector<int> vidx(nVert);
  int dim = 3;
  Eigen::VectorXd ue(24);
  int i = cellIdx / (N * N);
  cellIdx -= i*N*N;
  int j = cellIdx / N;
  cellIdx -= j*N;
  int k = cellIdx;
  for (int v = 0; v < nVert; v++){
    int ui = i + oneV[v][0];
    int uj = j + oneV[v][1];
    int uk = k + oneV[v][2];
    int uidx = ui*(N + 1)*(N + 1) + uj*(N + 1) + uk;
    vidx[v] = uidx;
  }
  return vidx;
}

///@param N size of cell. size of verts is N+1.
Eigen::VectorXd cell_u(int cellIdx, int N, const Eigen::VectorXd & u)
{
  int dim = 3;
  Eigen::VectorXd ue(24);
  std::vector<int> vidx = cell_verts(cellIdx, N);
  for (int v = 0; v < 8; v++){
    for (int d = 0; d < dim; d++){
      ue[v * dim + d] = u[3 * vidx[v] + d];
    }
  }
  return ue;
}

void testCopyU()
{
  //test copyPeriodicU
  Eigen::VectorXd u(81);
  Eigen::VectorXd u0(24);
  Eigen::Matrix3d strain = Eigen::Matrix3d::Identity();
  for (int i = 0; i < 24; i++){
    u0[i] = i;
  }
  copyPeriodicU(3, 1.0/3.0, u, u0,strain);
  for (int i = 0; i < u0.size(); i++){
    std::cout << i << " " << u0[i] << "\n";
  }
  std::cout << "========\n";
  for (int i = 0; i < u.size(); i++){
    std::cout << i << " " << u[i] << "\n";
  }
}

void test_cell_u()
{
  Eigen::VectorXd u(81);
  for (int i = 0; i < u.size(); i++){
    u[i] = i;
  }
  Eigen::VectorXd ue = cell_u(7, 2, u);
  for (int i = 0; i < ue.size(); i++){
    std::cout << i << " " << ue[i] << "\n";
  }
}

void FEM3DFun::setParam(const Eigen::VectorXd & x0)
{
  Vector3S color0(0.8f, 0.8f, 1.0f);
  param = x0;

  //read material distribution from field
  for (int ii = 0; ii < field->param.size(); ii++){
    field->setParam(ii, x0[ii]);
    //std::cout << x0[ii] << "\n";
  }
  
  for (int ii = 0; ii < gridSize[0]; ii++){
    for (int jj = 0; jj < gridSize[1]; jj++){
      for (int kk = 0; kk < gridSize[2]; kk++){
        Eigen::VectorXd coord = Eigen::VectorXd::Zero(3);
        int eIdx = gridToLinearIdx(ii, jj, kk, gridSize);
        coord[0] = (ii + 0.5) / gridSize[0];
        coord[1] = (jj + 0.5) / gridSize[1];
        coord[2] = (kk + 0.5) / gridSize[2];
        distribution[eIdx] = Clamp(field->f(coord), (T)1e-3, (T)1);
        (*fem.variable_coef)(ii, jj, kk) = femCoef(distribution[eIdx]);
      }
    }
  }
  fem.u.setZero(); fem.f.setZero(); fem.psi_D.clear(); fem.nodal_forces.clear();
  fem.Update_Matrices();
  hom->Homogenize();
  int n = 6;
  G = Eigen::MatrixXd::Zero(n, n);
  int i = 0;
  int j = 0;
  for (int i = 0; i < n; i++){
    copyPeriodicU(gridSize[0] + 1, fem.lattice.dx,u[i], hom->unit_test_u[i], hom->unit_test_strains[i]);
  }
  int nCell = gridSize[0] * gridSize[1] * gridSize[2];
  for (int cIdx = 0; cIdx < nCell; cIdx++){
    double coef = femCoef(distribution[cIdx]);
    for (int i = 0; i < n; i++){
      VectorX cell_u_i = cell_u(cIdx, gridSize[0], u[i]);
      //for (int j = 0; j < cell_u_i.size(); j++){
      //  std::cout << j << " " << cell_u_i[j] << "\n";
      //}
      //std::cout << "============\n";
      for (int j = i; j < n; j++){
        VectorX cell_u_j = cell_u(cIdx, gridSize[0], u[j]);
        int mat_id = 0;
        T e = cell_u_i.dot(Ke0[mat_id] * cell_u_j);
        G(i, j) += coef * e;
      }
    }
  }
  double vol = std::pow(lattice.dx * gridSize[0],3);
  G *= (1 / vol);
  //for (int j = 0; j < 6; j++){
  //  for (int k = 0; k < 6; k++){
  //    std::cout << G(j, k) << " ";
  //  }
  //  std::cout << "=========\n";
  //}
  //solve linear statics problems
  //density objective
  density = 0;
  for (unsigned int ii = 0; ii < distribution.size(); ii++){
    density += distribution[ii];
  }
  density /= distribution.size();

  Eigen::VectorXd X;
}

double FEM3DFun::f()
{
  double val = 0;
  Eigen::MatrixXd diff = G - G0;
  diff = diff.cwiseProduct(diff);
  double sum = 0;
  for (int i = 0; i < G.rows(); i++){
    for (int j = 0; j < G.cols(); j++){
      sum += wG(i, j) * diff(i, j) * diff(i, j);
    }
  }
  val = wG(0,0) * G(0, 0);// 0.5 * sum;
  val += 0.5 * mw * (density - m0) * (density - m0);
  std::cout << G(0,0)<<" "<<G(0,1)<<" "<<G(0,2) << "\n";
  std::cout << val << "\n";
  return val ;
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

  for (unsigned int ii = 0; ii < dfddist.size(); ii++){
    dfddist[ii] += (mw / distribution.size()) * (density - m0);
  }

  //volume of entire domain.
  double vol = std::pow(lattice.dx * gridSize[0], 3);
  int nEle = gridSize[0] * gridSize[1] * gridSize[2];
  int mat_id = 0;
  for (int row = 0; row < G.rows(); row++){
    for (int col = 0; col < G.cols(); col++){
      double dfdG = wG(row, col);// *(G(row, col) - G0(row, col)) / vol;
      for (int i = 0; i < nEle; i++){
        Eigen::VectorXd ui = cell_u(i, gridSize[0], u[row]);
        Eigen::VectorXd uj = cell_u(i, gridSize[0], u[col]);
        double dKdx = dfemCoef(distribution[i]);
        double dGdxi = dKdx * ui.dot(Ke0[mat_id] * uj);
        dfddist[i] += dfdG * dGdxi;
      }
    }
  }
  for (int i = 0; i < nEle; i++){
    std::cout<<dfddist[i]<<"\n";
  }
  Eigen::VectorXd coord = Eigen::VectorXd::Zero(3);
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
  out << G << "\n";
}

FEM3DFun::FEM3DFun() :m_Init(false),
gridSize(3,0),
field(0)
{
}

FEM3DFun::~FEM3DFun(){}
