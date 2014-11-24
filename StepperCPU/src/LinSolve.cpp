#include "LinSolve.hpp"
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"
#include "MatrixXd.hpp"
void linSolve(MatrixXd & AA, double * bb)
{
  int N = AA.mm;
  int info=0;
 // info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', N, N, 1, AA.M, N, bb, 1);
  if(info != 0){
    std::cout<<"Lapack solve error: "<<info<<"\n";
  }
}