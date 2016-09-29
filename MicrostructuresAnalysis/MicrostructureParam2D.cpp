#include "MicrostructureParam2D.hpp"

template <> void drawRectangle(const Rectangle2D & r, std::vector<double> & arr,
  const std::vector<int> & arrSize, double val);

template <int> void drawRectangle(const Rectangle2D & r, std::vector<int> & arr,
  const std::vector<int> & arrSize, int val);

template <>
void saveStructure2D(const std::vector<int> & arr,
  const std::vector<int> & arrSize, std::ostream & out);

template <>
void loadStructure2D(std::vector<double> & arr,
  std::vector<int> & arrSize, std::istream & in);

template <>
void make2DCubic(std::vector<double> & arr, const std::vector<int> & arrSize);

template <>
void make2DCubic(std::vector<int> & arr, const std::vector<int> & arrSize);

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


void test2DParam()
{
  int N = 32;
  int nEle = N * N;
  std::vector<double> arr(nEle, 1);
  std::vector<int> arrSize(2, N);
  
  //center rectangle
  Rectangle2D r;
  r.center[0] = 0.5;
  r.center[1] = 0.5;
  r.length[0] = 0.4;
  r.length[1] = 0.4;

  //branches on diagonal of center
  Rectangle2D r1;
  r1.center[0] = 0.25;
  r1.center[1] = 0.25;
  r1.length[0] = 0.2;
  r1.length[1] = 0.1;
  r1.rot = 3.1415/4;

  //soft beams on sides of center.
  Rectangle2D r2;
  r2.center[0] = 0.5;
  r2.center[1] = 0;
  r2.length[0] = 0.2;
  r2.length[1] = 0.5;

  //crosses
  Rectangle2D r3;
  r3.center[0] = 0;
  r3.center[1] = 0;
  r3.length[0] = 0.6;
  r3.length[1] = 0.2;

  //a test ellipse
  //Ellipse2D e;
  //e.center[0] = 0.5;
  //e.center[1] = 0.5;
  //e.r[0] = 0.4;
  //e.r[1] = 0.1;
  //e.rot = 3.14/3;
  //drawEllipse(e, arr, arrSize, 0.0);

  drawRectangle(r, arr, arrSize, 0.0);
  drawRectangle(r1, arr, arrSize, 0.0);
  drawRectangle(r2, arr, arrSize, 0.0);
  drawRectangle(r3, arr, arrSize, 0.0);

  make2DCubic(arr, arrSize);
  std::ofstream out("struct2d.txt");
  saveStructure2D(arr, arrSize, out);
  out.close();
}
