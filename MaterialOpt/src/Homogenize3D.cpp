#include "cfgDefs.h"
#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "Homogenize3D.hpp"
#include "pardiso_sym.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include "MaterialQuad.hpp"
#include "StrainLin.hpp"

typedef Eigen::Triplet<cfgScalar> TripletS;

void set_idx(std::vector<bool> & a, const std::vector<int> & idx, bool val)
{
  for (size_t ii = 0; ii < idx.size(); ii++){
    a[idx[ii]] = val;
  }
}

void saveArr(std::ostream & o, const std::vector<int> & a)
{
  o << a.size() << "\n";
  for (size_t i = 0; i < a.size(); i++){
    o << a[i] << "\n";
  }
}

Eigen::SparseMatrix<float>
makePeriodicSystem(const Eigen::SparseMatrix<float> & K0,
  const Eigen::SparseMatrix<float> & M43,
  const std::vector<int> & d1,
  const std::vector<int> & d2,
  const std::vector<int> & d3,
  const std::vector<int> & d4
  )
{
  Eigen::SparseMatrix<float> K22 = submatrix(K0, d2, d2, 1);
  Eigen::SparseMatrix<float> K23 = submatrix(K0, d2, d3, 1);
  Eigen::SparseMatrix<float> K24 = submatrix(K0, d2, d4, 1);
  Eigen::SparseMatrix<float> K33 = submatrix(K0, d3, d3, 1);
  Eigen::SparseMatrix<float> K34 = submatrix(K0, d3, d4, 1);
  Eigen::SparseMatrix<float> K44 = submatrix(K0, d4, d4, 1);

  Eigen::SparseMatrix<float> A12 = K23;
  A12 += K24*M43;
  Eigen::SparseMatrix<float> Ktop = concatRow(K22, A12);
  Eigen::SparseMatrix<float> A22 = K33;
  Eigen::SparseMatrix<float> KM = K34*M43;
  A22 += KM;
  A22 += Eigen::SparseMatrix<float>(KM.transpose());
  A22 += M43.transpose() * K44 * M43;
  Eigen::SparseMatrix<float> A21 = A12.transpose();
  Eigen::SparseMatrix<float> Kbottom = concatRow(A21, A22);
  Eigen::SparseMatrix<float> A = concatCol(Ktop, Kbottom);
  
  //debug output to matlatb.
  std::ofstream mout("A.txt");
  write_CSR(mout, A);
  mout.close();
  mout.open("K.txt");
  write_CSR(mout, K0);
  mout.close();
  mout.open("d.txt");
  mout << 4 << "\n";
  saveArr(mout, d1);
  saveArr(mout, d2);
  saveArr(mout, d3);
  saveArr(mout, d4);
  mout.close();
  mout.open("M43.txt");
  write_CSR(mout, M43);
  mout.close();
  mout.open("K22.txt");
  write_CSR(mout, K22);
  mout.close(); 
  
  return A;
  //return K22;
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
  
  d1 = expandIdx(n1, dim);
  
  //build M43 matrix that maps between free and fixed opposite faces and edge.
  int nEdgeVert = vertGridSize[0] + vertGridSize[1] + vertGridSize[2] - 6;
  int nEdgeDof = dim * nEdgeVert;
  int nFaceVert = ((vertGridSize[0] - 2) * (vertGridSize[1] - 2) + (vertGridSize[0] - 2) * (vertGridSize[2] - 2) +
    (vertGridSize[1] - 2) * (vertGridSize[2] - 2));
  int nFaceDof = dim * nFaceVert;
  int M43rows = 3 * nEdgeDof + nFaceDof;
  int M43cols = nEdgeDof + nFaceDof;
  M43 = Eigen::SparseMatrix<float>(M43rows, nEdgeDof + nFaceDof);
  M43.reserve(Eigen::VectorXi::Constant(M43cols, 3));
  for (int i = 0; i < nEdgeVert; i++){
    for (int j = 0; j < 3; j++){
      M43.insert(9 * i + j, 3 * i + j) = 1;
      M43.insert(9 * i + j + 3, 3 * i + j) = 1;
      M43.insert(9 * i + j + 6, 3 * i + j) = 1;
    }
  }
  for (int i = 0; i < nFaceDof; i++){
    M43.insert(i + 3 * nEdgeDof, i + nEdgeDof) = 1;
  }
 
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

  int nInternalVert = (vertGridSize[0] - 2) * (vertGridSize[1] - 2) * (vertGridSize[2] - 2);
  int nInternalDof = dim * nInternalVert;
  std::vector<int> n2(nInternalVert, 0);
  int cnt = 0;
  for (int i = 1; i < vertGridSize[0] - 1; i++){
    for (int j = 1; j < vertGridSize[1] - 1; j++){
      for (int k = 1; k < vertGridSize[2] - 1; k++){
        n2[cnt] = gridToLinearIdx(i, j, k, vertGridSize);
        cnt++;
      }
    }
  }
  d2 = expandIdx(n2, dim);

  //compute corner displacements U1
  //six boundary conditions.
  //number of equations to solve.
  int nEq = 6;
  Eigen::MatrixXf e0 = Eigen::MatrixXf::Identity(nEq, nEq);
  ufixed = Eigen::MatrixXf::Zero(d1.size(), e0.cols());
  float dx = em->X[1][2] - em->X[0][2];
  for (int j = 0; j < nEq; j++){
    Eigen::Matrix3f strain;
    strain << e0(j,0),     e0(j,5) / 2, e0(j,4) / 2,
              e0(j,5) / 2, e0(j,1),     e0(j,3) / 2,
              e0(j,4) / 2, e0(j,3) / 2, e0(j,2);
    ufixed.block(3, j, 3, 1) = dx * strain * Eigen::Vector3f(0.0f, 0.0f, (float)gridSize[2]);
    ufixed.block(6, j, 3, 1) = dx * strain * Eigen::Vector3f(0.0f, (float)gridSize[1], 0.0f);
    ufixed.block(9, j, 3, 1) = dx * strain * Eigen::Vector3f(0.0f, (float)gridSize[1], (float)gridSize[2]);
    ufixed.block(12, j, 3, 1) = dx * strain * Eigen::Vector3f((float)gridSize[0], 0.0f, 0.0f);
    ufixed.block(15, j, 3, 1) = dx * strain * Eigen::Vector3f((float)gridSize[0], 0.0f, (float)gridSize[2]);
    ufixed.block(18, j, 3, 1) = dx * strain * Eigen::Vector3f((float)gridSize[0], (float)gridSize[1], 0.0f);
    ufixed.block(21, j, 3, 1) = dx * strain * Eigen::Vector3f((float)gridSize[0], (float)gridSize[1], (float)gridSize[2]);
  }

  //Compute replative displacements of constrained vertices W = strain * y
  wfixed = Eigen::MatrixXf::Zero(d4.size(), nEq);
  //first loop over edges
  for (int i = 0; i < nEdgeVert; i++){
    //dof on the minus side
    int vm = n3[i];
    for (int j = 0; j < nEq; j++){
      Eigen::Matrix3f strain;
      strain << e0(j, 0), e0(j, 5) / 2, e0(j, 4) / 2,
        e0(j, 5) / 2, e0(j, 1), e0(j, 3) / 2,
        e0(j, 4) / 2, e0(j, 3) / 2, e0(j, 2);
      
      for (int k = 0; k < 3; k++){
        int n4idx = 3 * i + k;
        int vp = n4[n4idx];
        Eigen::Vector3f disp = em->X[vp] - em->X[vm];
        wfixed.block(n4idx * 3, j, 3, 1) = strain * disp;
      }
    }
  }

  //loop over faces
  for (int i = 0; i < nFaceVert; i++){
    int vm = n3[i + nEdgeVert];
    int n4idx = i + 3 * nEdgeVert;
    int vp = n4[n4idx];
    for (int j = 0; j < nEq; j++){
      Eigen::Matrix3f strain;
      strain << e0(j, 0), e0(j, 5) / 2, e0(j, 4) / 2,
        e0(j, 5) / 2, e0(j, 1), e0(j, 3) / 2,
        e0(j, 4) / 2, e0(j, 3) / 2, e0(j, 2);
      Eigen::Vector3f disp = em->X[vp] - em->X[vm];
      wfixed.block(n4idx * 3, j, 3, 1) = strain * disp;
    }
  }
  
  //Construct K22 K23+K24 K33+K34+K34'+K44
  bool triangular = false;
  bool fix_rigid = false;
  std::vector<StrainLin> ene(1);
  //E = 1e3.
  ene[0].param[0] = 100;
  ene[0].param[1] = 2400;
  std::vector<MaterialQuad * > material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    material[ii] = new MaterialQuad();
    for (unsigned int jj = 0; jj < material[ii]->e.size(); jj++){
      material[ii]->e[jj] = &ene[ii];
    }
    em->addMaterial(material[ii]);
  }
  em->initArrays();
  em->check();
  K0 = em->getStiffness(0);
  if (distribution.size() == 0){
    distribution = Eigen::VectorXd::Ones(em->e.size());
  }
  em->stiffnessPattern(m_I, m_J, triangular, fix_rigid);
  int rows = (int)m_I.size() - 1;
  std::vector<float> val(m_J.size(), 1);
  Eigen::SparseMatrix<float> pattern = Eigen::MappedSparseMatrix<float, Eigen::ColMajor>
    (rows, rows, m_J.size(), m_I.data(), m_J.data(), val.data());
  Eigen::SparseMatrix<float> K = makePeriodicSystem(pattern, M43, d1, d2, d3, d4);
  triangular = true;
  sparseToIJ(m_I, m_J, K, triangular);
  ////convert to 1-based
  checkSparseIndex(m_I, m_J);
  add(m_I, 1);
  add(m_J, 1);
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), K.rows(), pardisoState);
}

void Homogenize3D::solve()
{
  bool triangular = true;
  //number of rows in linear system
  int nrows = (int)m_I.size() - 1;
  //number of dof in the mesh;
  int ndof = dim * em->x.size();

  std::vector<int> vertGridSize = gridSize;
  for (size_t ii = 0; ii < vertGridSize.size(); ii++){
    vertGridSize[ii] ++;
  }
  int nEdgeVert = vertGridSize[0] + vertGridSize[1] + vertGridSize[2] - 6;
  int nEdgeDof = dim * nEdgeVert;
  int nFaceVert = ((vertGridSize[0] - 2) * (vertGridSize[1] - 2) + (vertGridSize[0] - 2) * (vertGridSize[2] - 2) +
    (vertGridSize[1] - 2) * (vertGridSize[2] - 2));
  int nFaceDof = dim * nFaceVert;
  //Compute full stiffness matrix
  getStiffnessSparse();
  //K for lhs
  Eigen::SparseMatrix<float> Kl = makePeriodicSystem(m_K, M43, d1, d2, d3, d4);
  sparseToVal(m_val, Kl, triangular);
  pardisoNumericalFactorize(m_I.data(), m_J.data(), m_val.data(), nrows, pardisoState);
  Eigen::SparseMatrix<float> K21 = submatrix(m_K, d2, d1, 1);
  Eigen::SparseMatrix<float> K31 = submatrix(m_K, d3, d1, 1);
  Eigen::SparseMatrix<float> K41 = submatrix(m_K, d4, d1, 1);
  Eigen::SparseMatrix<float> K24 = submatrix(m_K, d2, d4, 1);
  Eigen::SparseMatrix<float> K34 = submatrix(m_K, d3, d4, 1);
  Eigen::SparseMatrix<float> K44 = submatrix(m_K, d4, d4, 1);
  Eigen::SparseMatrix<float> Kbot = K31 + M43.transpose() * K41;
  Eigen::SparseMatrix<float> Umat = concatCol(K21, Kbot);
  Kbot = K34 + M43.transpose() * K44;
  Eigen::SparseMatrix<float> Wmat = concatCol(K24, Kbot);
  std::ofstream mout("W.txt");
  write_CSR(mout, Wmat);
  mout.close();
  mout.open("U.txt");
  write_CSR(mout, Umat);
  mout.close();
  //solve
  //|K22            K23 + K24M43                         | = - |  K21          | U1 - |   K24         | W
  //|(K23+K24M43)'  K33 + M34K43 + K43'M43 + M43'K44M43  |     | K31 + M43'K41 |      | K34 + M43'K44 |
  u.resize(ufixed.cols());
  Eigen::SparseMatrix<double> Kld = Kl.cast<double>();
  for (int i = 0; i < ufixed.cols(); i++){
    Eigen::VectorXf rhs = -Umat * ufixed.col(i) - Wmat * wfixed.col(i);
    if (i == 0){
      mout.open("rhs.txt");
      saveArr(rhs, mout);
      mout.close();
    }
    Eigen::VectorXd b = rhs.cast<double>();
    Eigen::VectorXd x = Eigen::VectorXd::Zero(nrows);
    pardisoBackSubstitute(m_I.data(), m_J.data(), m_val.data(), nrows, x.data(), b.data(), pardisoState);
    Eigen::VectorXd residual = Kld * x - b;
    double maxR = 0;
    for (int j = 0; j < residual.size(); j++){
      maxR = std::max(maxR, residual[j]);
    }
    std::cout << "max residual " << maxR << "\n";
    u[i].resize(ndof);
    //u(d1) = ufixed
    for (int j = 0; j < (int)d1.size(); j++){
      u[i][d1[j]] = ufixed(j,i);
    }
    //u(d2 d3) = x
    for (int j = 0; j < (int)d2.size(); j++){
      u[i][d2[j]] = x[j];
    }
    for (int j = 0; j < (int)d3.size(); j++){
      u[i][d3[j]] = x[j+d2.size()];
    }
    //u(d4) = U(d3) + wfixed.
    //edges followed by faces
    for (int j = 0; j < nEdgeVert; j++){
      for (int k = 0; k < 3; k++){
        for (int di = 0; di < 3; di++){
          int dof = 3 * (3 * j + k) + di;
          int vp = d4[dof];
          u[i][vp] = wfixed(dof, i) + x[d2.size() + 3*j+di];
        }
      }
    }
    for (int j = 0; j < nFaceDof; j++){
      int vp = d4[j + 3 * nEdgeDof];
      u[i][vp] = wfixed(3 * nEdgeDof + j, i) + x[j + nEdgeDof + d2.size()];
    }
  }
}

MatrixXS Homogenize3D::getKe(int ei)
{
  return (cfgScalar)(distribution[ei] / (3 - 2 * distribution[ei])) * 0.5* (K0+K0.transpose());
}

void Homogenize3D::getStiffnessSparse()
{
  int dim = 3;
  int N = dim * (int)em->x.size();
  std::vector<TripletS> coef;
  m_K=Eigen::SparseMatrix<cfgScalar> (N, N);
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
            cfgScalar val = K(dim * jj + dim1, dim * kk + dim2);
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
  m_K.setFromTriplets(coef.begin(), coef.end());
}