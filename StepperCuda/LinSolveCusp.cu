#include "LinSolveCusp.hpp"
#include "cuda.h"

#include <thrust/version.h>
#include <cusp/version.h>
#include <iostream>

#include <cusp/monitor.h>
#include <cusp/gallery/poisson.h>
#include <cusp/linear_operator.h>
#include <cusp/krylov/cg.h>

#include <cusp/coo_matrix.h>
#include <cusp/print.h>

#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>
#include <thrust/iterator/zip_iterator.h>
using namespace cusp;
LinSolveCusp::LinSolveCusp():device_I(0), device_J(0), device_V(0),
device_x(0), device_b(0), mSize(0), nnz(0), nIter(1000)
{}

void LinSolveCusp::init(std::vector<int> & I, std::vector<int> & J)
{
  mSize = I.size() - 1;
  nnz = I[I.size() - 1];
  cudaError_t status;
  status = cudaMalloc((void**)(&device_I), I.size() * sizeof(int));
  status = cudaMalloc((void**)(&device_J), nnz * sizeof(int));
  status = cudaMalloc((void**)(&device_V), nnz * sizeof(ValueType));
  status = cudaMalloc((void**)(&device_x), mSize * sizeof(ValueType));
  status = cudaMalloc((void**)(&device_b), mSize * sizeof(ValueType));

  status = cudaMemcpy(device_I, &(I[0]), I.size()*sizeof(int), cudaMemcpyHostToDevice);
  status = cudaMemcpy(device_J, &(J[0]), nnz*sizeof(int), cudaMemcpyHostToDevice);
}

void LinSolveCusp::solve(std::vector<ValueType> & A, ValueType * xIO)
{
  cudaMemcpy(device_V, &(A[0]), nnz*sizeof(ValueType), cudaMemcpyHostToDevice);
  cudaMemcpy(device_b, xIO, mSize*sizeof(ValueType), cudaMemcpyHostToDevice);
  std::vector <ValueType>zero(mSize, 0.0f);
  cudaMemcpy(device_x, &(zero[0]), mSize*sizeof(ValueType), cudaMemcpyHostToDevice);
  // wrap the device memory with a csr_matrix_view
  // *NOTE* raw pointers must be wrapped with thrust::device_ptr!
  thrust::device_ptr<int>   wrapped_device_I(device_I);
  thrust::device_ptr<int>   wrapped_device_J(device_J);
  thrust::device_ptr<ValueType> wrapped_device_V(device_V);
  thrust::device_ptr<ValueType> wrapped_device_x(device_x);
  thrust::device_ptr<ValueType> wrapped_device_b(device_b);
  
  {
    // use array1d_view to wrap the individual arrays
    typedef typename cusp::array1d_view< thrust::device_ptr<int>   > DeviceIndexArrayView;
    typedef typename cusp::array1d_view< thrust::device_ptr<ValueType> > DeviceValueArrayView;
    DeviceIndexArrayView row_offsets(wrapped_device_I, wrapped_device_I + mSize + 1);
    DeviceIndexArrayView column_indices(wrapped_device_J, wrapped_device_J + nnz);
    DeviceValueArrayView values(wrapped_device_V, wrapped_device_V + nnz);
    DeviceValueArrayView x(wrapped_device_x, wrapped_device_x + mSize);
    DeviceValueArrayView b(wrapped_device_b, wrapped_device_b + mSize);

    typedef cusp::csr_matrix_view < DeviceIndexArrayView,
      DeviceIndexArrayView,
      DeviceValueArrayView > DeviceView;
    DeviceView A(mSize, mSize, nnz, row_offsets, column_indices, values);
    cusp::default_monitor<ValueType> monitor(b, nIter, 1e-5);
    cusp::krylov::cg(A, x, b, monitor);
    cudaMemcpy((void*)xIO, device_x, mSize*sizeof(ValueType), cudaMemcpyDeviceToHost);
  }
}

LinSolveCusp::~LinSolveCusp()
{
  if (device_I != 0){
    cudaFree(device_I);
    device_I = 0;
  }

  if (device_J != 0){
    cudaFree(device_J);
    device_J = 0;
  }

  if (device_V != 0){
    cudaFree(device_V);
    device_V = 0;
  }

  if (device_x != 0){
    cudaFree(device_x);
    device_x = 0;
  }

  if (device_b != 0){
    cudaFree(device_b);
    device_b = 0;
  }
}

// where to perform the computation
typedef cusp::device_memory MemorySpace;
// which floating point type to use
typedef float ValueType;

void testCG(void)
{
  cusp::csr_matrix<int, ValueType, MemorySpace> P;
  // create a 2d Poisson problem on a 10x10 mesh
  cusp::gallery::poisson5pt(P, 100, 100);
  // allocate storage for solution (x) and right hand side (b)
  cusp::array1d<ValueType, MemorySpace> x(P.num_rows, 0);
  cusp::array1d<ValueType, MemorySpace> b(P.num_rows, 1);
  cusp::default_monitor<ValueType> monitor(b, 1000, 1e-3);
  // set preconditioner (identity)
  cusp::identity_operator<ValueType, MemorySpace> M(P.num_rows, P.num_rows);
  // solve the linear system A * x = b with the Conjugate Gradient method
  cusp::krylov::cg(P, x, b, monitor, M);

  // report solver results
  if (monitor.converged())
  {
    std::cout << "Solver converged to " << monitor.tolerance() << " tolerance";
    std::cout << " after " << monitor.iteration_count() << " iterations";
    std::cout << " (" << monitor.residual_norm() << " final residual)" << std::endl;
  }
  else
  {
    std::cout << "Solver reached iteration limit " << monitor.iteration_limit() << " before converging";
    std::cout << " to " << monitor.tolerance() << " tolerance ";
    std::cout << " (" << monitor.residual_norm() << " final residual)" << std::endl;
  }

  const int kb = 1024;
  const int mb = kb * kb;
  std::cout << "NBody.GPU" << std::endl << "=========" << std::endl << std::endl;

  std::cout << "CUDA version:   v" << CUDART_VERSION << std::endl;
  std::cout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << std::endl << std::endl;

  int devCount;
  cudaGetDeviceCount(&devCount);
  std::cout << "CUDA Devices: " << std::endl << std::endl;

  for (int i = 0; i < devCount; ++i)
  {
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, i);
    std::cout << i << ": " << props.name << ": " << props.major << "." << props.minor << std::endl;
    std::cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << std::endl;
    std::cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << std::endl;
    std::cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << std::endl;
    std::cout << "  Block registers: " << props.regsPerBlock << std::endl << std::endl;

    std::cout << "  Warp size:         " << props.warpSize << std::endl;
    std::cout << "  Threads per block: " << props.maxThreadsPerBlock << std::endl;
    std::cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]" << std::endl;
    std::cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]" << std::endl;
    std::cout << std::endl;
  }
  
}

void solveTriplet(std::vector<int> & I, std::vector<int> & J,
  std::vector<float> &val, float * x)
{
  int mSize = (int)I.size()-1;
  int nnz = I[I.size() - 1];
  cusp::csr_matrix<int, float, cusp::host_memory> A(mSize,mSize,nnz);
  cusp::array1d<int, cusp::host_memory> b(mSize);

  for (unsigned int ii = 0; ii < I.size(); ii++){
    A.row_offsets[ii] = I[ii];
  }
  for (unsigned int ii = 0; ii < J.size(); ii++){
    A.column_indices[ii] = J[ii];
    A.values[ii] = val[ii];
  }
  for (int ii = 0; ii < mSize; ii++){
    b[ii] = x[ii];
 //   std::cout << x[ii] << " " << "\n";
  }

  cusp::array1d<ValueType, MemorySpace> xdev(mSize, 0);
  cusp::array1d<ValueType, MemorySpace> bdev(b);
  cusp::csr_matrix<int, ValueType, MemorySpace> Adev(A);
  cusp::default_monitor<ValueType> monitor(bdev, 10000, 1e-7);
  cusp::identity_operator<ValueType, cusp::device_memory> M(A.num_rows, A.num_rows);
  cusp::krylov::cg(Adev, xdev, bdev, monitor, M);
 // cusp::print(A);
  std::cout << "converged: " << monitor.converged() << " " << monitor.residual_norm()<<"\n";
  cusp::array1d<ValueType, cusp::host_memory> result(xdev);
  //cusp::copy(xdev, b);
  //cusp::print(xdev);
  for (int ii = 0; ii < mSize; ii++){
    x[ii] = result[ii];
  }
  
}
