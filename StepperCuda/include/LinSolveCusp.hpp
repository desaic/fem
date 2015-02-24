#ifndef LIN_SOLVE_CUSP_HPP
#define LIN_SOLVE_CUSP_HPP
#include <vector>
class LinSolveCusp{
public:

};

void testCG();
void solveTriplet(std::vector<int> & I, std::vector<int> & J,
  std::vector<float> &val, float * x);
#endif