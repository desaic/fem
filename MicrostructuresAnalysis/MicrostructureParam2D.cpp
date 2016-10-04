#include "MicrostructureParam2D.hpp"
#include <fstream>

std::vector<int> linearToGridIdx(int linearIdx, const std::vector<int> & arrSize)
{
  int dim = (int)(arrSize.size());
  std::vector<int> coord(dim, 0);
  for (int i = 0; i < dim-1; i++){
    coord[i] = linearIdx / arrSize[i];
    linearIdx -= coord[i] * arrSize[i];
  }
  coord[dim - 1] = linearIdx;
  return coord;
}

int linearIdx(int x, int y, const std::vector<int> & arrSize)
{
  std::vector<int> coord(2);
  coord[0] = x;
  coord[1] = y;
  int idx = gridToLinearIdx(coord, arrSize);
  return idx;
}

int gridToLinearIdx(const std::vector<int> & coord, const std::vector<int> & arrSize)
{
  int idx = 0;
  int dim = (int)arrSize.size();
  int power = 1;
  for (int i = dim - 1; i >= 0; i--){
    idx += power * coord[i];
    power *= arrSize[i];
  }
  return idx;
}

bool inside(const std::vector<int> & coord,
  const std::vector<int> & gridSize)
{
  int dim = (int)gridSize.size(); 
  for (int i = 0; i < dim; i++){
    if (coord[i] < 0 || coord[i] >= gridSize[i]){
      return false;
    }
  }
  return true;
}

bool ptInRect(const float * pt, const float * len, int dim)
{
  for (int i = 0; i < 2; i++){
    if (pt[i] < - len[i] / 2
      || pt[i] > len[i] / 2){
      return false;
    }
  }
  return true;
}

bool ptInEllipse(const float * pt, const float * r, int dim)
{
  float sum = 0;
  sum = std::pow(pt[0] * r[1], 2.0f) + std::pow(pt[1] * r[0], 2.0f);
  float thresh = std::pow(r[0] * r[1], 2.0f);
  return sum <= thresh;
}
