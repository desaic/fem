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
  
  //find - and + face indices and interior indices d3 d4
  //in x y z order faces
  for (int ii = 0; ii < vertGridSize[1]; ii++){
    for (int jj = 0; jj < vertGridSize[2]; jj++){
      if ((ii == 0 || ii == vertGridSize[1] - 1) && (jj == 0 || jj == vertGridSize[2] - 1)){
        continue;
      }
      d3.push_back(3 * gridToLinearIdx(0, ii, jj, vertGridSize));
      d4.push_back(3 * gridToLinearIdx(vertGridSize[0] - 1, ii, jj, vertGridSize));
    }
  }
  for (int ii = 0; ii < vertGridSize[0]; ii++){
    for (int jj = 0; jj < vertGridSize[2]; jj++){
      if ((ii == 0 || ii == vertGridSize[0] - 1) && (jj == 0 || jj == vertGridSize[2] - 1)){
        continue;
      }
      d3.push_back(3 * gridToLinearIdx(ii, 0, jj, vertGridSize) + 1);
      d4.push_back(3 * gridToLinearIdx(ii, vertGridSize[1] - 1, jj, vertGridSize) + 1);
    }
  }
  for (int ii = 0; ii < vertGridSize[0]; ii++){
    for (int jj = 0; jj < vertGridSize[1]; jj++){
      if ((ii == 0 || ii == vertGridSize[0] - 1) && (jj == 0 || jj == vertGridSize[1] - 1)){
        continue;
      }
      d3.push_back(3 * gridToLinearIdx(ii, jj, 0, vertGridSize) + 2);
      d4.push_back(3 * gridToLinearIdx(ii, jj, vertGridSize[2] - 1, vertGridSize) + 2);
    }
  }
  
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

  //Compute positive face displacements W = strain * y

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
