#include "ElementMesh.hpp"
#include "Eigen/Sparse"
#include "MatrixXd.hpp"
void ElementMesh::getStiffnessSparse(std::vector<int> & I, std::vector<int> & J, 
                                     std::vector<float> &val)
{
  int N = 3* (int)x.size();
  Eigen::SparseMatrix<float> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    MatrixXd K  = getStiffness(ii);

  }
}
