#ifndef CONJUGATEGRADIENTCUDA_HPP
#define CONJUGATEGRADIENTCUDA_HPP
struct cublasContext;
struct cusparseContext;
struct cusparseMatDescr;
class ConjugateGradientCuda
{
public:
  ConjugateGradientCuda();
  ///@brief initialize gpu memories assuming matrix sparseness pattern is fixed
  int initCuda(int N,int _nz, int * _I, int * _J);
  int clearCuda();
  int solve( float * val, float * x, float * rhs);
  float tol;
  int max_iter;
  int devID;

//private:
  struct cublasContext* cublasHandle;
  struct cusparseContext * cusparseHandle;
  struct cusparseMatDescr * descr;
  
  int N,nz;
  int * I, * J;
  int *d_col, *d_row;
  float *d_val, *d_x;
  float *d_r, *d_p, *d_Ax;
};
#endif