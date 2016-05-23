#ifndef EIGENUTIL_HPP
#define EIGENUTIL_HPP

#include <Eigen/Sparse>
#include <vector>

void clampVector(Eigen::VectorXd & x, const Eigen::VectorXd & lb, const Eigen::VectorXd & ub);

///@brief write sparse matrix in matlab format such that one can copy paste
///into matlab workspace
void write_matlab(std::ostream &output, const char *variable_name,
                  const Eigen::SparseMatrix<double> & M);

///@brief write sparse matrix in matlab format for loading with matlab code.
void write_matlab_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M);

///@brief write sparse matrix in Vega format for loading with vega sparsematrix.
template <typename T>
void write_CSR(std::ostream &output,
              const Eigen::SparseMatrix<T> & M);

///@brief zero out off-diagonal terms of row and col ii if fixed[ii]!=0.
template <typename T>
void zeroOffDiag(Eigen::SparseMatrix<T> & K, const std::vector<int> & fixed)
{
  for(unsigned int col =0; col<K.cols(); col++){
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
      int row = it.row();
      if( (fixed[col] != 0 || fixed[row]!=0) && row != col){
        it.valueRef() = 0;
      }
    }
  }
}

void sparseToIJ(std::vector<int> & I,
  std::vector<int> & J, const Eigen::SparseMatrix<float> & K,
  bool triangular = false);

void sparseToVal(std::vector<double> & val,
  const Eigen::SparseMatrix<float> & K,
  bool triangular);

template <typename T>
Eigen::SparseMatrix<T>
rmConstrained(const Eigen::SparseMatrix<T> & K,
              const std::vector<int> & fixed)
{
  int cnt = 0;
  std::vector<int> cidx(K.rows(),0);
  for(unsigned int ii = 0; ii<fixed.size(); ii++){
    if(fixed[ii]){
      continue;
    }
    cidx[ii] = cnt;
    cnt++;
  }

  int N = cnt;
//  std::cout<<N<<"\n";
  Eigen::SparseMatrix<T> Kc(N, N);

  int maxnnz = 0;
  auto I = K.outerIndexPtr();
  for (auto ii = 1; ii <= K.cols(); ii++){
    maxnnz = std::max(I[ii] - I[ii-1], maxnnz);
  }
  maxnnz = std::min(maxnnz, nrow);
  Kc.reserve(Eigen::VectorXi::Constant(N, maxnnz));

  for(unsigned int col =0; col<K.cols(); col++){
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
      int row = it.row();
      if(fixed[col] || fixed[row]){
        continue;
      }
      Kc.insert(cidx[row], cidx[col]) = it.value();
    }
  }

  return Kc;
}

double sum(const Eigen::VectorXd & x);

std::vector<int>
expandIdx(const std::vector<int> & I,
int block_size);

///@param subrow subset of indices.
///@brief assuming column major matrix. Which is default of Eigen.
template <typename T>
Eigen::SparseMatrix<T>submatrix(const Eigen::SparseMatrix<T> & K,
  const std::vector<int> & subrow,
  const std::vector<int> & subcol,
  int block_size);

///@brief |A B|.
template <typename T>
Eigen::SparseMatrix<T>concatRow(const Eigen::SparseMatrix<T> & A,
  const Eigen::SparseMatrix<T> & B)
{
  Eigen::SparseMatrix<T> C(A.rows(), A.cols() + B.cols());
  if (A.rows() != B.rows()){
    std::cout << "Concat Row wrong rows " << A.rows() << " " << B.rows() << "\n";
    return C;
  }
  C.middleCols(0, A.cols()) = A;
  C.middleCols(A.cols(), B.cols()) = B;
  return C;
}

///@brief |A|
///       |B|.
template <typename T>
Eigen::SparseMatrix<T>concatCol(const Eigen::SparseMatrix<T> & A,
  const Eigen::SparseMatrix<T> & B)
{
  int nrow = A.rows() + B.rows();
  int maxnnz = 0;
  auto Ia = A.outerIndexPtr();
  auto Ib = B.outerIndexPtr();
  for (auto ii = 1; ii <= A.cols(); ii++){
    int nnz = Ia[ii] - Ia[ii - 1] + Ib[ii] - Ib[ii - 1];
    maxnnz = std::max(nnz, maxnnz);
  }
  maxnnz = std::min(maxnnz, nrow);
  
  Eigen::SparseMatrix<T> C(nrow, A.cols());
  C.reserve(Eigen::VectorXi::Constant(A.cols(), maxnnz));
  if (A.cols() != B.cols()){
    std::cout << "Concat Col wrong cols " << A.cols() << " " << B.cols() << "\n";
    return C;
  }
  
  for (int col = 0; col < A.cols(); col++){
    for (Eigen::SparseMatrix<T>::InnerIterator it(A, col); it; ++it){
      C.insert(it.row(), it.col())=it.value();
    }
    for (Eigen::SparseMatrix<T>::InnerIterator it(B, col); it; ++it){
      C.insert(it.row() + A.rows(), it.col())=it.value();
    }
  }

  return C;
}

#endif // EIGENUTIL_HPP
