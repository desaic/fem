#ifndef LIN_SOLVE_CUSP_HPP
#define LIN_SOLVE_CUSP_HPP
#include <vector>
class LinSolveCusp{
public:
  typedef float ValueType;
  
  LinSolveCusp();
  ~LinSolveCusp();
  void init(std::vector<int> & I, std::vector<int> & J);
  ///@param x Used as input and output.
  void solve(std::vector<ValueType> & A, ValueType * x);
  void dealloc();
  int * device_I, * device_J;
  ValueType * device_V, *device_x, *device_b;
  ///@brief matrix size assuming square matrix
  int mSize,nnz, nIter;
};

void testCG();
void solveTriplet(std::vector<int> & I, std::vector<int> & J,
  std::vector<float> &val, float * x);
#endif