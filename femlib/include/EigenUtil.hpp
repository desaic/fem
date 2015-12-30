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
void write_vega_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M);


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

  int max_nnz = 0;
  for(unsigned int col =0; col<K.cols(); col++){
    int n = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
      n++;
    }
    if(n>max_nnz){
      max_nnz = n;
    }
  }
  Kc.reserve(Eigen::VectorXi::Constant(N, max_nnz));

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

#endif // EIGENUTIL_HPP
