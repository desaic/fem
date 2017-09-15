
#include "NumericalCoarseningUsingMultigrid.h"

#include "MaterialQuad.hpp"
#include "StrainEne.hpp"

#include "cfgMaterialUtilities.h"

#include <iostream>

static const int d = 3;
Type_Define_VectorD_Outside_Class(d); Type_Define_VectorDi_Outside_Class(d); Type_Define_MatrixD_Outside_Class(d);

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

///@brief only works for cube. gridSize must be equal.
void swapDispXY(Eigen::VectorXd& u, int gridSize[3])
{
  int dim = 3;
  int N = gridSize[0]+1;
  int nVert = u.size() / dim;
  std::vector<int> vertSize(3, N);
  //swap displacements between vertices.
  for (int i = 0; i < N-1; i++){
    for (int j = i+1; j < N; j++){
      for (int k = 0; k < N; k++){
        int v0Idx = cfgMaterialUtilities::getGridToVectorIndex(i, j, k, vertSize[0], vertSize[1], vertSize[2]);
        int v1Idx = cfgMaterialUtilities::getGridToVectorIndex(j, i, k, vertSize[0], vertSize[1], vertSize[2]);
        for (int l = 0; l < dim; l++){
          double tmp = u[dim * v0Idx + l];
          u[dim * v0Idx + l] = u[dim * v1Idx + l];
          u[dim * v1Idx + l] = tmp;
        }
      }
    }
  }
  //swap x y coordinates for each displacement
  for (int i = 0; i < nVert; i++){
    double tmp = u[i*dim];
    u[i*dim] = u[i*dim+1];
    u[i*dim+1] = tmp;
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





NumericalCoarseningUsingMultigrid::NumericalCoarseningUsingMultigrid()
{
  m_nbBlockRep = 1;
  m_nbSubdiv = 1;

  m_dim = d;
  m_hom = NULL;

  m_useVariableDensity = true;
  m_densityMin = 1.e-4;

  m_init = false;
}

NumericalCoarseningUsingMultigrid::~NumericalCoarseningUsingMultigrid()
{
  SAFE_DELETE(m_hom);
}

void NumericalCoarseningUsingMultigrid::init(int N[3], int iNbBlockRep, int iNbSubdiv, const std::vector<MaterialQuad> &iBaseMaterials, StructureType iType)
{
  m_init = true;

  m_dim = 3;
  m_nbBlockRep = iNbBlockRep;
  m_nbSubdiv = iNbSubdiv;

  int n[3] = {N[0]*m_nbBlockRep*m_nbSubdiv, N[1]*m_nbBlockRep*m_nbSubdiv, N[2]*m_nbBlockRep*m_nbSubdiv};

  m_baseMaterials3D = iBaseMaterials;

  std::vector<int> cellMaterials(n[0]*n[1]*n[2], 0);

  const int solver_type = 1;
  const int processor_type = 1;
  const bool write_residual = false;
  const int num_mg_level = 5;
  const int des_iter_num = 1;
  const int asc_iter_num = 1;
  const int v_cycle_iter_num = 1;
  const int krylov_iter_num = 1000;
  int material_pattern = 0;

  T length = (T)1; VectorDi cell_counts = VectorDi::Ones()*(n[0] - 1);
  Lattice<d> lattice(cell_counts, length / (T)(cell_counts[0] + 1));
  Lattice<d> lattice_output = Lattice<d>(VectorDi::Ones()*n[0], length / (T)(cell_counts[0] + 1));
  m_fem.Initialize(lattice);
  m_fem.solver_type = ElasticHexFEM<d>::SolverType(solver_type);
  m_fem.processor_type = ElasticHexFEM<d>::ProcessorType(processor_type);
  m_fem.des_iter_num = des_iter_num;
  m_fem.asc_iter_num = asc_iter_num;
  m_fem.v_cycle_iter_num = v_cycle_iter_num;
  m_fem.krylov_iter_num = krylov_iter_num;
  m_fem.use_implicit_psi_P = true;
  m_fem.num_multigrid_level = num_mg_level;

  ////reinitialize materials
  int nmat = (int)iBaseMaterials.size();
  m_fem.materials.clear();
  for (int imat=0; imat<nmat; imat++)
  {
    cfgScalar mu = iBaseMaterials[imat].e[0]->param[0];
    cfgScalar lambda = iBaseMaterials[imat].e[0]->param[1];

    cfgScalar E=0, nu=0;
    cfgMaterialUtilities::fromLamesParametersToYoungModulusPoissonRatio(lambda, mu, E, nu);

    m_fem.materials.push_back(ElasticMaterial((T)E, (T)nu));
  }
  //if (m_useVariableDensity)
  {
    std::swap(m_fem.materials[0], m_fem.materials[1]);
  }

  for (LatticeIterator<d> iter(lattice, LatticeIterator<d>::CELL); iter.Valid(); iter.Next())
  {
    m_fem.material_id(iter.Coord()) = 0;
  }

  m_fem.u.setZero(); m_fem.f.setZero(); m_fem.psi_D.clear(); m_fem.nodal_forces.clear();
  m_u.clear();

  m_fem.Update_Matrices();
  m_K0 = m_fem.Ke;

  int ndisp = 6;
  m_u.resize(6);
  int nVert = (n[0]+1)*(n[1]+1)*(n[2]+1);
  for (int i = 0; i < ndisp; i++)
  {
    m_u[i].resize(d*nVert);
  }

  Eigen::Vector3i varSize( n[0],  n[1], n[2]);
  SAFE_DELETE(m_hom);
  m_hom = new Homogenization<d>(m_fem);
  m_hom->type = (iType==Cubic? Homogenization<d>::CUBIC: Homogenization<d>::DEFAULT);
  m_hom->Initialize();
  if (m_useVariableDensity)
  {
    m_fem.variable_coef = new Field<T, d>();
    m_fem.variable_coef->Resize(varSize); 
    m_fem.variable_coef->Fill((T)0);
  }
}

void NumericalCoarseningUsingMultigrid::computeCoarsenedElasticityTensorAndParameters3D(const std::vector<int> &iMaterials, int N[3], const std::vector<cfgScalar> &iBaseMaterialDensities, StructureType iType, 
                                                                                        std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters)
{
  std::vector<int> cellMaterials;
  cfgMaterialUtilities::repMaterialAssignment(N[0], N[1], N[2], iMaterials, m_nbBlockRep, m_nbBlockRep, m_nbBlockRep, m_nbSubdiv, cellMaterials); 

  int n[3] = {N[0]*m_nbBlockRep*m_nbSubdiv, N[1]*m_nbBlockRep*m_nbSubdiv, N[2]*m_nbBlockRep*m_nbSubdiv};

  if (m_useVariableDensity)
  {
    for (int i=0; i<n[0]; i++)
    {
      for (int j=0; j<n[1]; j++)
      {
        for (int k=0; k<n[2]; k++)
        {
          int indMat = cfgMaterialUtilities::getGridToVectorIndex(i, j, k, n[0], n[1], n[2]);
          int mat_id = cellMaterials[indMat];
          (*m_fem.variable_coef)(i, j, k) = (mat_id==1? 1: m_densityMin);
        }
      }
    }
  }
  else
  {
    for (int i=0; i<n[0]-1; i++)
    {
      for (int j=0; j<n[1]-1; j++)
      {
        for (int k=0; k<n[2]-1; k++)
        {
          int indMat = cfgMaterialUtilities::getGridToVectorIndex(i, j, k, n[0], n[1], n[2]);
          int mat_id = cellMaterials[indMat];
          /*if (mat_id==1)
          {
            mat_id = 0;
          }
          else
          {
            mat_id = -1; //void
          }*/ 
          m_fem.material_id(i, j, k) = mat_id;
        }
      }
    }
  }
 
  m_fem.u.setZero(); m_fem.f.setZero(); m_fem.psi_D.clear(); m_fem.nodal_forces.clear();
  assert(m_hom);
  std::streambuf* orig_buf = std::cout.rdbuf();
  std::cout.rdbuf(NULL);
  m_hom->Homogenize();
  std::cout.rdbuf(orig_buf);

  int ndisp = 6;
  if (m_hom->type == Homogenization<d>::CUBIC)
  {
    copyPeriodicU(n[0] + 1, m_fem.lattice.dx, m_u[0], m_hom->unit_test_u[0], m_hom->unit_test_strains[0]);
    m_u[1] = m_u[0];
    swapDispXY(m_u[1], n);
    copyPeriodicU(n[0] + 1, m_fem.lattice.dx, m_u[3], m_hom->unit_test_u[1], m_hom->unit_test_strains[1]);
  }
  else{
    for (int i = 0; i < ndisp; i++)
    {
      copyPeriodicU(n[0] + 1, m_fem.lattice.dx, m_u[i], m_hom->unit_test_u[i], m_hom->unit_test_strains[i]);
    }
  }

  Eigen::MatrixXd  C = Eigen::MatrixXd::Zero(ndisp, ndisp);
  int nCell = n[0] * n[1] * n[2];
  for (int cIdx = 0; cIdx < nCell; cIdx++)
  {
    for (int i = 0; i < ndisp; i++)
    {
      VectorX cell_u_i = cell_u(cIdx, n[0], m_u[i]);
      for (int j = i; j < ndisp; j++)
      {
        VectorX cell_u_j = cell_u(cIdx, n[0], m_u[j]);
        int ii, jj, kk;
        cfgMaterialUtilities::getVectorIndexToGrid(cIdx, n[0]-1, n[1]-1, n[2]-1, ii, jj, kk);
        int mat_id = m_fem.material_id(ii%(n[0]-1), jj%(n[1]-1), kk%(n[2]-1));
        if (mat_id>=0)
        {
          T e = cell_u_i.dot(m_K0[mat_id] * cell_u_j);
          cfgScalar coeff = (m_useVariableDensity? (cellMaterials[cIdx]==1? 1: m_densityMin): 1);
          C(i, j) += coeff*e;
        }
      }
    }
  }
  double vol = std::pow(m_fem.lattice.dx*n[0], 3);
  C *= (1 / vol);

  std::vector<float> Cvalues;
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
    {
      if (j<=i)
      {
        ioTensorValues.push_back(C(i,j));
      }
    }
  }

  if (iType == Cubic)
  {
    cfgScalar r = cfgMaterialUtilities::computeMaterialDensity(iMaterials, iBaseMaterialDensities);
    ioParameters.push_back(r);

    cfgScalar G = C(3, 3);
    cfgScalar a = C(0, 1) / C(0, 0);
    cfgScalar nu_xy = a / (1 + a);
    cfgScalar E = C(0, 1)*(1 + nu_xy)*(1 - 2 * nu_xy) / nu_xy;

    ioParameters.push_back(E);
    ioParameters.push_back(nu_xy);
    ioParameters.push_back(G);
  }
  else
  {
    assert(0); //not implemented
  }
}
