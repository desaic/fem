#include <iostream>
#include "EigenUtil.hpp"
#include "runtime.hpp"

void testEigenUtil()
{
  int N = 36;
  Eigen::SparseMatrix<double> M(N, N);
  M.reserve(Eigen::VectorXi::Constant(N, N / 2 + 1));
  int cnt = 1;
  for (int ii = 0; ii < N; ii++){
    int j0 = (ii % 2);
    for (int jj = j0; jj < N; jj += 2){
      M.insert(jj, ii) = cnt;
      cnt++;
    }
  }
  std::cout << M << "\n";
  int block_size = 3;
  int subsize = 4;
  std::vector<int> n1(subsize), n2(subsize);
  for (size_t ii = 0; ii < n1.size(); ii++){
    n1[ii] = std::min((int)ii * 4, N / block_size - 1);
    n2[ii] = std::min((int)ii * 4 + 1, N / block_size - 1);
  }

  Eigen::SparseMatrix<double> sub = submatrix(M, n1, n2, block_size);
  std::cout << sub << "\n";

  N = 3;
  cnt = 1;
  Eigen::SparseMatrix<double> A(N, N), B(2 * N, N), D(N, 2 * N);
  A.reserve(Eigen::VectorXi::Constant(N, N / 2 + 1));
  B.reserve(Eigen::VectorXi::Constant(N, N));
  B.reserve(Eigen::VectorXi::Constant(2 * N, N));
  for (int ii = 0; ii < N; ii++){
    int j0 = (ii % 2);
    for (int jj = j0; jj < N; jj += 2){
      A.insert(jj, ii) = cnt;
      cnt++;
    }
    for (int jj = j0; jj < 2 * N; jj += 2){
      B.insert(jj, ii) = cnt;
      cnt++;
    }
  }
  for (int ii = 0; ii < 2 * N; ii++){
    int j0 = (ii % 2);
    for (int jj = j0; jj < N; jj += 2){
      D.insert(jj, ii) = cnt;
      cnt++;
    }
  }
  std::cout << A << "\n" << B << "\n" << D << "\n";
  Eigen::SparseMatrix<double> C = concatRow(A, D);
  std::cout << C << "\n";
  C = concatCol(B, A);
  std::cout << C << "\n";
}