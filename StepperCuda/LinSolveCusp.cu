#include "LinSolveCusp.hpp"
#include "cuda.h"
#include <thrust/version.h>
#include <cusp/version.h>
#include <iostream>

#include <cusp/hyb_matrix.h>
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


  // allocate storage for (4,3) matrix with 6 nonzeros
  cusp::csr_matrix<int, float, cusp::host_memory> A(4, 3, 6);

  // initialize matrix entries on host
  A.row_offsets[0] = 0;  // first offset is always zero
  A.row_offsets[1] = 2;
  A.row_offsets[2] = 2;
  A.row_offsets[3] = 3;
  A.row_offsets[4] = 6; // last offset is always num_entries

  A.column_indices[0] = 0; A.values[0] = 10;
  A.column_indices[1] = 2; A.values[1] = 20;
  A.column_indices[2] = 2; A.values[2] = 30;
  A.column_indices[3] = 0; A.values[3] = 40;
  A.column_indices[4] = 1; A.values[4] = 50;
  A.column_indices[5] = 2; A.values[5] = 60;

  // A now represents the following matrix
  //    [10  0 20]
  //    [ 0  0  0]
  //    [ 0  0 30]
  //    [40 50 60]

  // print matrix entries
  cusp::print(A);


}
