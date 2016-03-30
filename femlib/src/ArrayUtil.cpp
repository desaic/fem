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
  std::cout << "nnz " << I[I.size() - 1] << " "<< J.size() << "\n";
  for (int ii = 0; ii < I.size() - 1; ii++){
    //pardiso needs diagonal even if 0
    bool hasDiagonal = false;
    for (int jj = I[ii]; jj < I[ii + 1]; jj++){
      if (J[jj] >= nrow){
        std::cout << ii << " " << jj << " "<<J[jj]<<  "\n";
        return -1;
      }
      if (J[jj] == ii){
        hasDiagonal = true;
      }
    }
    if (!hasDiagonal){
      std::cout << "Col " << ii << " has no diagonal.\n";
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

void add(std::vector<int> & a, int v)
{
  for (size_t i = 0; i < a.size(); i++){
    a[i] += v;
  }
}

void saveArr(const Eigen::VectorXf & v, std::ostream & out)
{
  out << v.size() << "\n";
  for (int i = 0; i < v.size(); i++){
    out << v[i] << "\n";
  }
}
