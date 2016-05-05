#ifndef ARRAYUTIL_HPP
#define ARRAYUTIL_HPP
#include <vector>
#include <Eigen/Dense>

#include <iostream>

void BBox(const std::vector<Eigen::Vector3f >& v,
  Eigen::Vector3f & mn, Eigen::Vector3f & mx);

template<typename T>
void add(std::vector<T> & dst, const std::vector<T> & src)
{
  for(unsigned int ii = 0;ii<dst.size();ii++){
    dst[ii] += src[ii];
  }
}

template<typename T>
std::vector<T> mul(float f, const std::vector<T> & src)
{
  std::vector<T> prod(src.size());
  for(unsigned int ii = 0;ii<src.size();ii++){
    prod[ii] = f*src[ii];
  }
  return prod;
}

template<typename T>
void addmul(std::vector<T> & dst, float f, const std::vector<T> & src)
{
  for(unsigned int ii = 0;ii<src.size();ii++){
    dst[ii] += f*src[ii];
  }
}

double infNorm(const Eigen::VectorXd & a);

int checkSparseIndex(const std::vector<int > & I, const std::vector<int> & J);

int gridToLinearIdx(int ix, int iy, int iz, const std::vector<int> & gridSize);

int gridToLinearIdx(int ix, int iy, const std::vector<int> & gridSize);

int linearIdx(const Eigen::VectorXi & idx, const std::vector<int> & gridSize);


///@brief add f to each 3 subvector of a.
void addVector3d(std::vector<double> & a, const Eigen::Vector3d & f,
  const std::vector<int> & idx);

void copyVert3(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<double> & u);

void copyVert3(Eigen::VectorXd & x, const std::vector<int> & vidx,
  const std::vector<Eigen::Vector3f> & X);
void add(std::vector<int> & a, int v);

void saveArr(const Eigen::VectorXf & v, std::ostream & out);
#endif // ARRAYUTIL_HPP
