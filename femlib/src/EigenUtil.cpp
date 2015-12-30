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

void write_vega_lists(std::ostream &output,
                        const Eigen::SparseMatrix<double> & M)
{
  output<<M.rows()<<"\n"<<M.cols()<<"\n";
  output.precision(16);
  for(int i=0; i<M.cols(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
      output<< i << " " << it.row() << " " << it.value() <<"\n";
    }
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
