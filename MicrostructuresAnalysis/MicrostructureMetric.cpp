#include "MicrostructureMetric.hpp"
#include "MicrostructureParam2D.hpp"
#include <cmath>

std::vector<double> gaussianFilter(int ksize, double sigma)
{
  std::vector<double> filter(ksize);
  int krad = ksize / 2;
  double sum = 0;
  for (int i = 0; i < ksize; i++){
    double x = i - krad;
    double G = std::exp(-0.5 * std::pow(x / sigma, 2));
    sum += G;
    filter[i] = G;
  }
  for (int i = 0; i < ksize; i++){
    filter[i] /= sum;
  }
  return filter;
}

double
MicrostructureMetric2D::d(const std::vector<double> & s1, 
  const std::vector<double> & s2,
  const std::vector<int> & gridSize)
{
  int ksize = 5;
  double sigma = 1;
  if (s1.size() != s2.size()){
    std::cout << "2D metric Error. s1 and s2 different size.\n";
    return -1;
  }
  std::vector<double> filter = gaussianFilter(ksize, sigma);
  std::vector<double> blur1 = s1;
  std::vector<double> blur2 = s2;

  blur2D(blur1, gridSize, filter);
  blur2D(blur2, gridSize, filter);
  double dist = 0;
  for (size_t i = 0; i < blur1.size(); i++){
    double diff = blur1[i] - blur2[i];
    dist += std::pow(diff, 2);
  }
  return dist;
}

void blur2D(std::vector<double> & s, const std::vector<int> gridSize,
  const std::vector<double> & filter)
{
  int filterRadius = (int)filter.size() / 2;
  //blur x
  std::vector<double> buf(gridSize[0]);
  for (int j = 0; j < gridSize[1]; j++){
    for (int i = 0; i < gridSize[0]; i++){
      int l = linearIdx(i, j, gridSize);
      buf[i] = s[l];
    }
    for (int i = 0; i < gridSize[0]; i++){
      double sum = 0;
      for (int k = 0; k < (int)filter.size(); k++){
        int xIdx = i - filterRadius + k;
        xIdx = std::max(0, xIdx);
        xIdx = std::min(gridSize[0] - 1, xIdx);
        sum += buf[xIdx] * filter[k];
      }
      int l = linearIdx(i, j, gridSize);
      s[l] = sum;
    }
  }

  //blur y
  buf.resize(gridSize[1]);
  for (int i = 0; i < gridSize[0]; i++){
    for (int j = 0; j < gridSize[1]; j++){  
      int l = linearIdx(i, j, gridSize);
      buf[j] = s[l];
    }
    for (int j = 0; j < gridSize[1]; j++){
      double sum = 0;
      for (int k = 0; k < (int)filter.size(); k++){
        int yIdx = j - filterRadius + k;
        yIdx = std::max(0, yIdx);
        yIdx = std::min(gridSize[1] - 1, yIdx);
        sum += buf[yIdx] * filter[k];
      }
      int l = linearIdx(i, j, gridSize);
      s[l] = sum;
    }
  }
}
