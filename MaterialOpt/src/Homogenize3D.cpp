#include "ArrayUtil.hpp"
#include "Homogenize3D.hpp"
#include "pardiso_sym.hpp"
#include <iostream>
#include <Eigen/Sparse>

template <typename T>
Eigen::SparseMatrix<T>submatrix(const Eigen::SparseMatrix<T> & K,
  const std::vector<int> & i1,
  const std::vector<int> & i2,
  int block_size)
{
  Eigen::SparseMatrix<T> sub(i1.size()*block_size, i2.size()*block_size);

  return sub;
}

void Homogenize3D::init()
{
  if (em == 0){
    std::cout << "Homogenize3D uninit mesh.\n";
    return;
  }

  pardisoState = new PardisoState();
  pardisoInit(pardisoState);

  //find corner vertices n1
  n1.resize(8, 0);
  n1[0] = gridToLinearIdx(0, 0, 0, gridSize);
  n1[1] = gridToLinearIdx(0, 0, gridSize[2]-1, gridSize);
  n1[2] = gridToLinearIdx(0, gridSize[1] - 1, 0, gridSize);
  n1[3] = gridToLinearIdx(0, gridSize[1] - 1, gridSize[2] - 1, gridSize);
  n1[4] = gridToLinearIdx(gridSize[0] - 1, 0, 0, gridSize);
  n1[5] = gridToLinearIdx(gridSize[0] - 1, 0, gridSize[2] - 1, gridSize);
  n1[6] = gridToLinearIdx(gridSize[0] - 1, gridSize[1] - 1, 0, gridSize);
  n1[7] = gridToLinearIdx(gridSize[0] - 1, gridSize[1] - 1, gridSize[2], gridSize);

  //find - and + face indices and interior indices n3 n4
    
  //compute corner displacements U1

  //Compute positive face displacements W = strain * y

  //Construct K22 K23+K24 K33+K34+K34'+K44


}

void Homogenize3D::solve()
{
  //Compute full stiffness matrix

  
  //solve
  //|K22            K23+K24        | = - |  K21    | U1 - |   K24   | W
  //|(K23+K24)'     K33+K34+K43+K44|     | K31+K41 |      | K34+K44 |
}
