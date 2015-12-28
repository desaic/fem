#include "ArrayUtil.hpp"
#include <iostream> 
void BBox(const std::vector<Eigen::Vector3f >& v,
  Eigen::Vector3f & mn, Eigen::Vector3f & mx)
{
  mn = v[0];
  mx = v[0];
  for(unsigned int ii = 1 ;ii<v.size();ii++){
    for(int dim = 0 ; dim<3;dim++){
      if(v[ii][dim]<mn[dim]){
        mn[dim] = v[ii][dim];
      }
      if(v[ii][dim]>mx[dim]){
        mx[dim] = v[ii][dim];
      }
    }
  }
}

int checkSparseIndex(const std::vector<int > & I, const std::vector<int> & J)
{
  int nrow = (int)I.size() - 1;
  int maxIdx = 0;
  for (int ii = 0; ii < I.size() - 1; ii++){
    for (int jj = I[ii]; jj < I[ii + 1]; jj++){
      maxIdx = std::max(J[ii], maxIdx);
      if (J[jj] >= nrow){
        std::cout << ii << " " << jj << "\n";
        return -1;
      }
    }
  }
  std::cout << "max idx " << maxIdx << "\n";
  return 0;
}
