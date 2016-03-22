#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "ElementMesh.hpp"
#include "Homogenize3D.hpp"
#include "pardiso_sym.hpp"
#include <iostream>
#include <Eigen/Sparse>

void set_idx(std::vector<bool> & a, const std::vector<int> & idx, bool val)
{
  for (size_t ii = 0; ii < idx.size(); ii++){
    a[idx[ii]] = val;
  }
}

void Homogenize3D::init()
{
  if (em == 0){
    std::cout << "Homogenize3D uninit mesh.\n";
    return;
  }
  
  pardisoState = new PardisoState();
  pardisoInit(pardisoState);

  //find groups of constrained and free vertices.
  //find corner vertices n1
  std::vector<int> n1(8, 0);
  std::vector<int> vertGridSize = gridSize;
  for (size_t ii = 0; ii < vertGridSize.size(); ii++){
    vertGridSize[ii] ++;
  }
  n1[0] = gridToLinearIdx(0, 0, 0, vertGridSize);
  n1[1] = gridToLinearIdx(0, 0, vertGridSize[2] - 1, vertGridSize);
  n1[2] = gridToLinearIdx(0, vertGridSize[1] - 1, 0, vertGridSize);
  n1[3] = gridToLinearIdx(0, vertGridSize[1] - 1, vertGridSize[2] - 1, vertGridSize);
  n1[4] = gridToLinearIdx(vertGridSize[0] - 1, 0, 0, vertGridSize);
  n1[5] = gridToLinearIdx(vertGridSize[0] - 1, 0, vertGridSize[2] - 1, vertGridSize);
  n1[6] = gridToLinearIdx(vertGridSize[0] - 1, vertGridSize[1] - 1, 0, vertGridSize);
  n1[7] = gridToLinearIdx(vertGridSize[0] - 1, vertGridSize[1] - 1, vertGridSize[2]-1, vertGridSize);
  
  std::vector<int> d1 = expandIdx(n1, dim);
  
  //build M43 matrix
  int nEdgeDof = dim * (vertGridSize[0] + vertGridSize[1] + vertGridSize[2] - 6);
  int nFaceDof = dim * ((vertGridSize[0] - 2) * (vertGridSize[1] - 2) + (vertGridSize[0] - 2) * (vertGridSize[2] - 2) +
    (vertGridSize[1] - 2) * (vertGridSize[2] - 2));
  int M43rows = 3 * nEdgeDof + nFaceDof;
  int M43cols = nEdgeDof + nFaceDof;
  M43 = Eigen::SparseMatrix<float>(M43rows, nEdgeDof + nFaceDof);
  M43.reserve(Eigen::VectorXi::Constant(M43cols, 3));
  for (int i = 0; i < nEdgeDof; i++){
    for (int j = 0; j < 3; j++){
      M43.insert(3 * i + j, i) = 1;
    }
  }
  for (int i = 0; i < nFaceDof; i++){
    M43.insert(i + 3 * nEdgeDof, i + nEdgeDof) = 1;
  }
  std::cout << M43 << "\n";
  
  //x y z edges
  //
  std::vector<int> n3, n4;
  for (int i = 1; i < vertGridSize[0] - 1; i++){
    n3.push_back(gridToLinearIdx(i, 0, 0, vertGridSize));
    n4.push_back(gridToLinearIdx(i, vertGridSize[1] - 1, 0, vertGridSize));
    n4.push_back(gridToLinearIdx(i, vertGridSize[1] - 1, vertGridSize[2] - 1, vertGridSize));
    n4.push_back(gridToLinearIdx(i, 0, vertGridSize[2] - 1, vertGridSize));
  }
  for (int i = 1; i < vertGridSize[1] - 1; i++){
    n3.push_back(gridToLinearIdx(0, i, 0, vertGridSize));
    n4.push_back(gridToLinearIdx(0, i, vertGridSize[2] - 1, vertGridSize));
    n4.push_back(gridToLinearIdx(vertGridSize[0] - 1, i, vertGridSize[2] - 1, vertGridSize));
    n4.push_back(gridToLinearIdx(vertGridSize[0] - 1, i, 0, vertGridSize));
  }
  for (int i = 1; i < vertGridSize[2] - 1; i++){
    n3.push_back(gridToLinearIdx(0, 0, i, vertGridSize));
    n4.push_back(gridToLinearIdx(vertGridSize[0] - 1, 0, i, vertGridSize));
    n4.push_back(gridToLinearIdx(vertGridSize[0] - 1, vertGridSize[1] - 1, i, vertGridSize));
    n4.push_back(gridToLinearIdx(0, vertGridSize[1] - 1, i, vertGridSize));
  }
  //find - and + face indices and interior indices d3 d4
  //in x y z order faces
  for (int ii = 1; ii < vertGridSize[1] - 1; ii++){
    for (int jj = 1; jj < vertGridSize[2] - 1; jj++){
      n3.push_back(gridToLinearIdx(0, ii, jj, vertGridSize));
      n4.push_back(gridToLinearIdx(vertGridSize[0] - 1, ii, jj, vertGridSize));
    }
  }
  for (int ii = 1; ii < vertGridSize[0] - 1; ii++){
    for (int jj = 1; jj < vertGridSize[2] - 1; jj++){
      n3.push_back(gridToLinearIdx(ii, 0, jj, vertGridSize));
      n4.push_back(gridToLinearIdx(ii, vertGridSize[1] - 1, jj, vertGridSize));
    }
  }
  for (int ii = 1; ii < vertGridSize[0] - 1; ii++){
    for (int jj = 1; jj < vertGridSize[1] - 1; jj++){
      n3.push_back(gridToLinearIdx(ii, jj, 0, vertGridSize));
      n4.push_back(gridToLinearIdx(ii, jj, vertGridSize[2] - 1, vertGridSize));
    }
  }
  d3 = expandIdx(n3, dim);
  d4 = expandIdx(n4, dim);

  int nVert = vertGridSize[0] * vertGridSize[1] * vertGridSize[2];
  int nDof = dim * nVert;
  std::vector<bool> usedDOF(nDof, false);
  set_idx(usedDOF, d1, true);
  set_idx(usedDOF, d3, true);
  set_idx(usedDOF, d4, true);
  for (int ii = 0; ii < nDof; ii++){
    if (!usedDOF[ii]){
      d2.push_back(ii);
    }
  }

  //compute corner displacements U1
  //six boundary conditions.
  Eigen::MatrixXf e0 = Eigen::MatrixXf::Identity(6, 6);
  ufixed = Eigen::MatrixXf::Zero(d1.size(), e0.cols());
  float dx = em->X[1][2] - em->X[0][2];
  for (int j = 0; j < e0.cols(); j++){
    Eigen::Matrix3f strain;
    strain << e0(j,0),     e0(j,5) / 2, e0(j,4) / 2,
              e0(j,5) / 2, e0(j,1),     e0(j,3) / 2,
              e0(j,4) / 2, e0(j,3) / 2, e0(j,2);
    ufixed.block(3, j, 3, 1) = strain * dx * Eigen::Vector3f(0.0f, 0.0f, (float)gridSize[2]);
    ufixed.block(6, j, 3, 1) = strain * dx * Eigen::Vector3f(0.0f, (float)gridSize[1], 0.0f);
    ufixed.block(9, j, 3, 1) = strain * dx * Eigen::Vector3f(0.0f, (float)gridSize[1], (float)gridSize[2]);
    ufixed.block(12, j, 3, 1) = strain * dx * Eigen::Vector3f((float)gridSize[0], 0.0f, 0.0f);
    ufixed.block(15, j, 3, 1) = strain * dx * Eigen::Vector3f((float)gridSize[0], 0.0f, (float)gridSize[2]);
    ufixed.block(18, j, 3, 1) = strain * dx * Eigen::Vector3f((float)gridSize[0], (float)gridSize[1], 0.0f);
    ufixed.block(21, j, 3, 1) = strain * dx * Eigen::Vector3f((float)gridSize[0], (float)gridSize[1], (float)gridSize[2]);
  }

  //Compute positive face displacements W = strain * y
  wfixed = Eigen::VectorXf::Zero(d4.size(), ufixed.cols());
  for (size_t ii = 0; ii < d4.size(); ii++){
    //dof on the plus side
    int dofp = d4[ii];
    int vertp = dofp / dim;
    int dofm = d3[ii];
    int vertm = dofm / dim;
    Eigen::Vector3f y = em->X[vertp] - em->X[vertm];
    for (int j = 0; j < ufixed.cols(); j++){
      
    }
  }
  //Construct K22 K23+K24 K33+K34+K34'+K44
  bool triangular = true;
  bool fix_rigid = false;
  em->stiffnessPattern(m_I, m_J, triangular, fix_rigid);
  int rows = (int)m_I.size() - 1;
  std::vector<double> val(m_J.size(), 1);
  Eigen::SparseMatrix<double> P = Eigen::MappedSparseMatrix<double, Eigen::ColMajor>
    (rows, rows, m_J.size(), m_I.data(), m_J.data(), val.data());
  Eigen::SparseMatrix<double> sub = submatrix(P, n1, n1, 3);
  std::cout << sub << "\n";
}

void Homogenize3D::solve()
{
  //Compute full stiffness matrix

  
  //solve
  //|K22            K23+K24        | = - |  K21    | U1 - |   K24   | W
  //|(K23+K24)'     K33+K34+K43+K44|     | K31+K41 |      | K34+K44 |
}
