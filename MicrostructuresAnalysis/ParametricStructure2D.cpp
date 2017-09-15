#include "ParametricStructure2D.hpp"
#include "MicrostructureMetric.hpp"
#include "FileUtil.hpp"
#include <iostream>

ParametricStructure2D::ParametricStructure2D(int Nx, int Ny)
{
  resize(Nx, Ny);
}

void
ParametricStructure2D::resize(int Nx, int Ny)
{
  if (Nx <= 0 || Ny <= 0){
    std::cout << "Parametric structure N<=0.\n";
    return;
  }
  gridSize.resize(2);
  gridSize[0] = Nx;
  gridSize[1] = Ny;
  int nEle = Nx * Ny;
  s.resize(nEle);
}

ParametricNPR::ParametricNPR()
{
  rec_center.center[0] = 0.5f;
  rec_center.center[1] = 0.5f;
  rec_center.length[0] = 0.4f;
  rec_center.length[1] = 0.4f;

  //branches on diagonal of center
  rec_branch.center[0] = 0.3f;
  rec_branch.center[1] = 0.3f;
  rec_branch.length[0] = 0.4f;
  rec_branch.length[1] = 0.1f;
  rec_branch.rot = 3.1415f / 4.0f;

  //soft beams on sides of center.
  rec_beam.center[0] = 0.5f;
  rec_beam.center[1] = 0.0f;
  rec_beam.length[0] = 0.1f;
  rec_beam.length[1] = 0.4f;

  //crosses
  rec_cross.center[0] = 0.0f;
  rec_cross.center[1] = 0.0f;
  rec_cross.length[0] = 0.75f;
  rec_cross.length[1] = 0.15f;

  rec_cross_tip.center[0] = 0.3f;
  rec_cross_tip.center[1] = 0.0f;
  rec_cross_tip.r[0] = 0.05f;
  rec_cross_tip.r[1] = 0.15f;
  rec_cross_tip.rot = 0;
  extractParam();
}

//distribute parameters into geometric primitives.
void ParametricNPR::distributeParam()
{
  rec_center.length[0] = (float)param[0];
  rec_center.length[1] = (float)param[0];

  rec_branch.length[0] = (float)param[1];
  rec_branch.length[1] = (float)param[2];

  rec_beam.length[0] = (float)param[3];
  rec_beam.length[1] = (float)param[4];

  rec_cross.length[0] = (float)param[5];
  rec_cross.length[1] = (float)param[6];
}

//read variable parameters from geometric primitives
void ParametricNPR::extractParam()
{
  param.resize(7);
  param[0] = rec_center.length[0];
  
  param[1] = rec_branch.length[0];
  param[2] = rec_branch.length[1];

  param[3] = rec_beam.length[0];
  param[4] = rec_beam.length[1];

  param[5] = rec_cross.length[0];
  param[6] = rec_cross.length[1];
}

// center rectangle size with equal x and y size.
// branch length
// branch width
// beam length
// beam width
// cross length
// cross width
void
ParametricNPR::eval(const std::vector<double> & p)
{
  param = p;
  distributeParam();

  std::fill(s.begin(), s.end(), 1);
  drawRectangle(rec_center, s, gridSize, 0.0);
  drawRectangle(rec_branch, s, gridSize, 0.0);
  drawRectangle(rec_beam, s, gridSize, 0.0);
  drawRectangle(rec_cross, s, gridSize, 0.0);
  //drawEllipse(rec_cross_tip, s, gridSize, 0.0);
  make2DCubic(s, gridSize);
}


template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

//using coordinate descent.
void fitTarget(ParametricStructure2D * st , const std::vector<double> & target,
  const std::vector<int> & gridSize)
{
  int nParam = (int)st->param.size();
  float step0 = 1e-2;
  MicrostructureMetric2D metric;
  int iter = 0;
  int maxIter = 10;
  std::vector<double> p = st->param;
  double eps = 1e-4;
  int nSearch = 10;
  
  for (int iter = 0; iter < maxIter; iter++){
    bool change = false;
    for (int i = 0; i < nParam; i++){
      std::vector<double> p0 = p;
      double f0 = metric.d(st->s, target, gridSize);
      bool found = false;
      float h = 0;
      double fp = 0, fm;
      for (int j = 0; j < nSearch; j++){
        // find minimum step size that changes the objective.
        h = step0 * std::pow(2,j);
        p[i] = p0[i] + h;
        st->eval(p);
        fp = metric.d(st->s, target, gridSize);
        if (std::abs(fp - f0) > eps){
          found = true;
          break;
        }
      }
      
      if (!found){
        //this param converged.
        continue;
      }
      p[i] = p0[i] - h;
      st->eval(p);
      fm = metric.d(st->s, target, gridSize);
      if ( (f0 <= fm+eps) && (f0 <= fp+eps)) {
        //current point is a local minimum.
        continue;
      }
      change = true;
      double grad = (fp - fm) / (2 * h);
      //line search
      double sign = sgn(grad);
      double h0 = h;
      double fmin = f0;
      double hmin = 0;
      for (int j = 0; j < nSearch; j++){
        h = - sign * h0 * std::pow(2, j);
        p[i] = p0[i] + h;
        st->eval(p);
        fp = metric.d(st->s, target, gridSize);
        if (fp > fmin){
          break;
        }
        else{
          fmin = fp;
          hmin = h;
        }
      }
      p[i] = p0[i] + hmin;
    }
    if (!change){
      break; 
    }
  }
}

void test2DParam()
{
  int N = 128;
  std::vector<int> gridSize(2, N);
  ParametricNPR npr;
  npr.resize(N, N);
  npr.eval(npr.param);
  std::ofstream out("struct2d.txt");
  std::vector<double> filter = gaussianFilter(5, 1);
  std::vector<double> target;
  std::vector<int> inputSize;
  std::ifstream in("target.txt");
  loadStructure2D(target, inputSize, in);
  std::vector<double> targetResize;
  resize2D(target, inputSize, targetResize, gridSize);
  fitTarget(&npr, targetResize, gridSize);
  saveStructure2D(npr.s, npr.gridSize, out);
  out.close();
}
