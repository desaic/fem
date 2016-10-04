#ifndef MICROSTRUCTURE_METRIC_HPP
#define MICROSTRUCTURE_METRIC_HPP

#include <vector>

class MicrostructureMetric2D{
public:
  double d(const std::vector<double> & s1, const std::vector<double> & s2, 
    const std::vector<int> & gridSize);
};

std::vector<double> gaussianFilter(int ksize, double sigma);

void blur2D(std::vector<double> & s, const std::vector<int> gridSize,
  const std::vector<double> & filter);

#endif