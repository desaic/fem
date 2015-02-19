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
  cusp::hyb_matrix<int, ValueType, MemorySpace> P;
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


  // dimensions of the matrix
  int num_rows = 3;
  int num_cols = 3;

  // number of (i,j,v) triplets
  int num_triplets = 10;

  // allocate storage for unordered triplets
  cusp::array1d<int, cusp::device_memory> I(num_triplets);  // row indices
  cusp::array1d<int, cusp::device_memory> J(num_triplets);  // column indices
  cusp::array1d<float, cusp::device_memory> V(num_triplets);  // values

  // fill triplet arrays
  I[0] = 2; J[0] = 0; V[0] = 10;
  I[1] = 0; J[1] = 2; V[1] = 10;
  I[2] = 1; J[2] = 1; V[2] = 10;
  I[3] = 2; J[3] = 0; V[3] = 10;
  I[4] = 1; J[4] = 1; V[4] = 10;
  I[5] = 0; J[5] = 0; V[5] = 10;
  I[6] = 2; J[6] = 2; V[6] = 10;
  I[7] = 0; J[7] = 0; V[7] = 10;
  I[8] = 1; J[8] = 0; V[8] = 10;
  I[9] = 0; J[9] = 0; V[9] = 10;

  // sort triplets by (i,j) index using two stable sorts (first by J, then by I)
  thrust::stable_sort_by_key(J.begin(), J.end(), thrust::make_zip_iterator(thrust::make_tuple(I.begin(), V.begin())));
  thrust::stable_sort_by_key(I.begin(), I.end(), thrust::make_zip_iterator(thrust::make_tuple(J.begin(), V.begin())));

  // compute unique number of nonzeros in the output
  int num_entries = thrust::inner_product(thrust::make_zip_iterator(thrust::make_tuple(I.begin(), J.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(I.end(), J.end())) - 1,
    thrust::make_zip_iterator(thrust::make_tuple(I.begin(), J.begin())) + 1,
    int(0),
    thrust::plus<int>(),
    thrust::not_equal_to< thrust::tuple<int, int> >()) + 1;

  // allocate output matrix
  cusp::coo_matrix<int, float, cusp::device_memory> A(num_rows, num_cols, num_entries);

  // sum values with the same (i,j) index
  thrust::reduce_by_key(thrust::make_zip_iterator(thrust::make_tuple(I.begin(), J.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(I.end(), J.end())),
    V.begin(),
    thrust::make_zip_iterator(thrust::make_tuple(A.row_indices.begin(), A.column_indices.begin())),
    A.values.begin(),
    thrust::equal_to< thrust::tuple<int, int> >(),
    thrust::plus<float>());

  // print matrix
  cusp::print(A);

}
