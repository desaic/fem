#include "EigenUtil.hpp"
#include <assert.h>
void clampVector(Eigen::VectorXd & x, const Eigen::VectorXd & lb, const Eigen::VectorXd & ub)
{
  if (lb.size() > 0){
    assert(lb.size() == x.size());
    for (int ii = 0; ii < x.size(); ii++){
      x[ii] = std::max(lb[ii], x[ii]);
    }
  }
  if (ub.size()>0){
    assert(ub.size() == x.size());
    for (int ii = 0; ii < x.size(); ii++){
      x[ii] = std::min(ub[ii], x[ii]);
    }
  }
}

void write_matlab(std::ostream &output, const char *variable_name,
                  const Eigen::SparseMatrix<double> & M)
{

  output<<"I=[";
   for(int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<i+1;
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }

   output<<"];\n  J=[";
   for(int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<it.row()+1;
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }

   output<<"];\n  V=[";
   for(int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<it.value();
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }
   output<<"];\n";
   output<<variable_name<<"=sparse(I,J,V,";
   output<<M.cols()<<", "<<M.rows()<<");"<<std::endl;
}

void write_matlab_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M)
{
  int nnz = M.nonZeros();
  output<<M.rows()<<" "<<nnz<<"\n";
  for(int i=0; i<M.cols(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
      output<< i+1 << " " << it.row()+1 << " " << it.value() <<"\n";
    }
  }
}

template <typename T>
void write_vega_lists(std::ostream &output,
  const Eigen::SparseMatrix<T> & M)
{
  output<<M.rows()<<"\n"<<M.cols()<<"\n";
  output.precision(12);
  for(int i=0; i<M.cols(); ++i){
    for (Eigen::SparseMatrix<T>::InnerIterator it(M, i); it; ++it){
      output<< i << " " << it.row() << " " << it.value() <<"\n";
    }
  }
}

void sparseToIJ(std::vector<int> & I,
  std::vector<int> & J, const Eigen::SparseMatrix<float> & K,
  bool triangular)
{
  I.resize(K.cols() + 1);
  J.clear();
  I[0] = 0;
  for (int col = 0; col < K.cols(); col++){
    for (Eigen::SparseMatrix<float>::InnerIterator it(K, col); it; ++it){
      if (triangular && it.row() > it.col()){
        continue;
      }
      J.push_back(it.row());
    }
    I[col + 1] = J.size();
  }
}

double sum(const Eigen::VectorXd & x)
{
  double sum = 0;
  for (int ii = 0; ii < x.size(); ii++){
    sum += x[ii];
  }
  return sum;
}

std::vector<int>
expandIdx(const std::vector<int> & I,
  int block_size)
{
  std::vector<int> J(I.size() * block_size);
  for (size_t ii = 0; ii < I.size(); ii++){
    for (int bb = 0; bb < block_size; bb++){
      J[ii*block_size + bb] = I[ii] * block_size + bb;
    }
  }
  return J;
}

///@param subrow subset of indices.
///@brief assuming column major matrix. Which is default of Eigen.
template <typename T>
Eigen::SparseMatrix<T>submatrix(const Eigen::SparseMatrix<T> & K,
  const std::vector<int> & _subrow,
  const std::vector<int> & _subcol,
  int block_size)
{
  std::vector<int> subrow = expandIdx(_subrow, block_size);
  std::vector<int> subcol = expandIdx(_subcol, block_size);
  int ncol = (int)subcol.size();
  int nrow = (int)subrow.size();
  Eigen::SparseMatrix<T> sub(nrow, ncol);
  std::vector<int> rowmap(K.rows(), -1);

  for (auto ii = 0; ii < subrow.size(); ii++){
    rowmap[subrow[ii]] = ii;
  }
  int maxnnz = 0;
  auto I = K.outerIndexPtr();
  for (auto ii = 0; ii < subcol.size(); ii++){
    int col = subcol[ii];
    maxnnz = std::max(I[col + 1] - I[col], maxnnz);
  }
  maxnnz = std::min(maxnnz, nrow);
  sub.reserve(Eigen::VectorXi::Constant(ncol, maxnnz));
  for (int cc = 0; cc < (int)subcol.size(); cc++){
    int col = subcol[cc];
    for (Eigen::SparseMatrix<T>::InnerIterator it(K, col); it; ++it){
      int row = it.row();
      int rowInSubmatrix = rowmap[row];
      if (rowInSubmatrix < 0){
        continue;
      }
      T val = it.value();
      sub.insert(rowInSubmatrix, cc) = val;
    }
  }
  return sub;
}

template 
Eigen::SparseMatrix<double>submatrix<double>(const Eigen::SparseMatrix<double> & K,
  const std::vector<int> & subrow,
  const std::vector<int> & subcol,
  int block_size);

template 
Eigen::SparseMatrix<float>submatrix<float>(const Eigen::SparseMatrix<float> & K,
  const std::vector<int> & subrow,
  const std::vector<int> & subcol,
  int block_size);

template
void write_vega_lists(std::ostream &output,
const Eigen::SparseMatrix<double> & M);

template
void write_vega_lists(std::ostream &output,
  const Eigen::SparseMatrix<float> & M);
