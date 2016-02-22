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

int gridToLinearIdx(int ix, int iy, int iz, const std::vector<int> & gridSize)
{
  return ix * gridSize[1] * gridSize[2] + iy * gridSize[2] + iz;
}

int gridToLinearIdx(int ix, int iy, const std::vector<int> & gridSize)
{
  return ix * gridSize[1] + iy;
}

int linearIdx(const Eigen::VectorXi & idx, const std::vector<int> & gridSize)
{
  return idx[0] * gridSize[1] * gridSize[2] + idx[1] * gridSize[2] + idx[2];
}



void copyVert3(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<Eigen::Vector3f> & X)
{
  int dim = 3;
  x = Eigen::VectorXd(dim * vidx.size());
  for (unsigned int ii = 0; ii < vidx.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      x[ii*dim + jj] = X[vidx[ii]][jj];
    }
  }
}

void copyVert3(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<double> & u)
{
  int dim = 3;
  x = Eigen::VectorXd(dim * vidx.size());
  for (unsigned int ii = 0; ii < vidx.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      x[ii*dim + jj] = u[dim*vidx[ii] + jj];
    }
  }
}

void addVector3d(std::vector<double> & a, const Eigen::Vector3d & f,
  const std::vector<int> & idx)
{
  int dim = (int)f.size();
  for (unsigned int ii = 0; ii < idx.size(); ii++){
    for (int jj = 0; jj < dim; jj++){
      a[dim * idx[ii] + jj] += f[jj];
    }
  }
}
