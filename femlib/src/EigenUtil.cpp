#include "EigenUtil.hpp"

void write_matlab(std::ostream &output, const char *variable_name,
                  const Eigen::SparseMatrix<double> & M)
{

  output<<"I=[";
   for(unsigned int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<i+1;
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }

   output<<"];\n  J=[";
   for(unsigned int i=0; i<M.cols(); ++i){
     for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
       output<<" "<<it.row()+1;
     }
     if(i>0 && i % 100==0){
       output<<" ...\n";
     }
   }

   output<<"];\n  V=[";
   for(unsigned int i=0; i<M.cols(); ++i){
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
  for(unsigned int i=0; i<M.cols(); ++i){
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
  for(unsigned int i=0; i<M.cols(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
      output<< i << " " << it.row() << " " << it.value() <<"\n";
    }
  }
}
