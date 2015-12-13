#include "ArrayUtil.hpp"
#include "ConfigFile.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "EigenUtil.hpp"
#include "Element.hpp"
#include "FileUtil.hpp"
#include "GridField2D.hpp"

#include "LinPardiso.hpp"

#include "StepperGrad.hpp"
#include "StepperStatic.hpp"
#include "StepperNewtonDyn.hpp"
#include "StepperNewmark.hpp"
#include "StepperLinNewmark.hpp"
#include "StepperReplay.hpp"

#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "StrainCorotLin.hpp"
#include "Quadrature.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <random>

typedef Eigen::Triplet<double> Tripletd;
typedef std::vector<std::vector<int > >  Grid2D;

std::vector<float> trajectory;

LinPardiso<double> linsolver;

double expo = 1.1;

double weight_M = 1e-6;
double Mfrac = 0.3;

double dx0 = 0.0045;
double dy0 = 0.0045;
double xweight = 5;
double yweight = 1;

Eigen::SparseMatrix<double> Kglobal;

void makeQuadGrid(ElementMesh * em, int nx, int ny, float meshScale = 1.0f);

std::vector<int> topVerts(ElementMesh * em, const Grid2D & grid);
std::vector<int> botVerts(ElementMesh * em, const Grid2D & grid);
std::vector<int> leftVerts(ElementMesh * em, const Grid2D & grid);
std::vector<int> rightVerts(ElementMesh * em, const Grid2D & grid);

std::vector<double> elementDisp(int eidx, ElementMesh * em);

int max_abs(std::vector<double> & a);

void loadArr2Di(std::istream & in, std::vector<std::vector<int> > &arr);

void clampVec(std::vector<double> &v, double mn, double mx);

void discreteSearch(ElementMesh * em, Stepper *stepper, const std::vector<std::vector<int> > & grid);

void fixRigid2d(Eigen::SparseMatrix<double> & K, ElementMesh * mesh);

void enforcePeriodicity(Eigen::SparseMatrix<double> & K, ElementMesh * mesh,
                        const std::vector<std::vector<int> > & grid);


double measureStretchY(ElementMesh * em, const std::vector<std::vector<int> > & grid)
{
  std::vector<int> tv, bv;
  tv = topVerts(em, grid);
  bv = botVerts(em, grid);
  double stretch = 0;
  for(unsigned int ii = 0; ii<tv.size(); ii++){
    stretch += (em->x[tv[ii]][1] - em->x[bv[ii]][1]) - (em->X[tv[ii]][1] - em->X[bv[ii]][1]);
  }
  stretch/=bv.size();
  return stretch;
}

double measureStretchX(ElementMesh * em, const std::vector<std::vector<int> > & grid)
{
  std::vector<int> lv, rv;
  lv = leftVerts(em, grid);
  rv = rightVerts(em, grid);
  double stretch = 0;
  for(unsigned int ii = 0; ii<lv.size(); ii++){
    stretch += (em->x[rv[ii]][0] - em->x[lv[ii]][0]) - (em->X[rv[ii]][0] - em->X[lv[ii]][0]);
  }
  stretch/=lv.size();
  return stretch;
}

double sum(const std::vector<double> & v)
{
  double sum = 0;
  for(unsigned int ii = 0; ii<v.size(); ii++){
    sum += v[ii];
  }
  return sum;
}

void setEntries(std::vector<double> & v, const std::vector<int>& fixedDof, double val)
{
  for(unsigned int jj = 0; jj<v.size(); jj++){
    if(fixedDof[jj]){
      v[jj] = val;
    }
  }
}

//poisson's ratio objective
double G(ElementMesh * em, const Grid2D & grid)
{
  double dx, dy;
  dx = measureStretchX(em, grid);
  dy = measureStretchY(em, grid);
  trajectory.push_back(dx);
  trajectory.push_back(dy);
  double val = 0.5 * xweight * (dx-dx0) * (dx-dx0) + 0.5 * yweight * (dy-dy0)*(dy-dy0);
  double massweight = 10;
  double mass = 0;
  for(int ii = 0; ii<em->depth.size(); ii++){
    mass += em->depth[ii];
  }
//  val += mass*massweight;
  return val;
}

std::vector<double> DGDu(ElementMesh * em, const std::vector<std::vector<int> > & grid)
{
  int N = em->X.size();
  int dim = em->dim;
  std::vector<double>grad(dim * N, 0);
  double dx = measureStretchX(em, grid);
  double dy = measureStretchY(em, grid);
  std::vector<int> verts;
  verts = topVerts(em, grid);
  for(unsigned int ii = 0; ii<verts.size(); ii++){
    grad[dim * verts[ii]+1] += (dy-dy0)*yweight;
  }
  verts = botVerts(em, grid);
  for(unsigned int ii = 0; ii<verts.size(); ii++){
    grad[dim * verts[ii]+1] -= (dy-dy0)*yweight;
  }
  verts = rightVerts(em, grid);
  for(unsigned int ii = 0; ii<verts.size(); ii++){
    grad[dim * verts[ii]] += (dx-dx0)*xweight;
  }
  verts = leftVerts(em, grid);
  for(unsigned int ii = 0; ii<verts.size(); ii++){
    grad[dim * verts[ii]] -= (dx-dx0)*xweight;
  }
  for(unsigned int ii = 0; ii<grad.size(); ii++){
    grad[ii] /=verts.size();
  }
  return grad;
}

void testdF(ElementMesh * em, const Eigen::VectorXd & Ke)
{
  std::vector<double> anaGrad = evaldF(em, Ke);
}

void addTimes(std::vector<double> & a, double c, const std::vector<double> & b){
  for(unsigned int ii = 0; ii<a.size(); ii++){
    a[ii] += c*b[ii];
  }
}

int max_abs(std::vector<double> & a)
{
  double maxval= 0;
  int maxidx = -1;
  for(unsigned int ii=0; ii<a.size(); ii++){
    double val = std::abs(a[ii]);
    if(val>maxval){
      maxval = val;
      maxidx = ii;
    }
  }
  return maxidx;
}

void scaleVec(std::vector<double>&v, double c)
{
  for(unsigned int ii = 0; ii<v.size(); ii++){
    v[ii] *= c;
  }
}

///@brief displacement vector of a single element
std::vector<double> elementDisp(int eidx, ElementMesh * em)
{
  int dim = em->dim;
  Element * ele = em->e[eidx];
  std::vector<double> u(dim * ele->nV());
  for(int ii = 0; ii<ele->nV(); ii++){
    int vidx = ele->at(ii);
    Eigen::Vector3d ui = em->x[vidx] - em->X[vidx];
    for(int jj = 0; jj<dim; jj++){
      u[dim * ii + jj] = ui[jj];
    }
  }
  return u;
}

void clampVec(std::vector<double> &v, double mn, double mx)
{
  for(unsigned int ii = 0; ii<v.size(); ii++){
    if(v[ii] > mx){
      v[ii] = mx;
    }
    if(v[ii]<mn){
      v[ii] =  mn;
    }
  }
}

void readGrid2D(GridField2D & matgrid, std::vector<double> & a, int nx, int ny)
{
  for(int ii = 0; ii<nx; ii++){
    for(int jj =0; jj<ny; jj++){
      Eigen::VectorXd x(2);
      x(0) = (ii+0.5)/nx;
      x(1) = (jj+0.5)/ny;
      double val = matgrid.f(x);
      a[ii*ny + jj]= val;
    }
  }
}

std::vector<double> DGDp(GridField2D & matgrid, const std::vector<double> & DGDx,
                 int nx, int ny)
{
  std::vector<double>grad(matgrid.param.size(), 0);
  for(int ii = 0; ii<nx; ii++){
    for(int jj =0; jj<ny; jj++){
      Eigen::VectorXd x(2);
      x(0) = (ii+0.5)/nx;
      x(1) = (jj+0.5)/ny;
      Eigen::SparseVector<double> DxDp = matgrid.df(x);
      int eidx = ii * ny + jj;
      for (Eigen::SparseVector<double>::InnerIterator it(DxDp); it; ++it){
        grad[it.index()] += DGDx[eidx] * it.value();
      }
    }
  }
  return grad;
}

void OptMat(World * world)
{
  Stepper * stepper = world->stepper;
  std::ofstream out_eval;
  out_eval.open("eval.txt");
  std::ofstream log("log.txt");
  //convert discrete material to continuous density
  double minDensity = 1e-3;
  double maxDensity = 1;
  int ITER = 0;
  double maxStep = 0.5;

  FileUtilIn in("../microstructure_SA.txt");
  std::vector<std::vector< int> > microstructures;
  loadArr2Di(in.in, microstructures);

  ElementMesh * em = world->em_[0];
  std::vector<std::vector<int> > grid;
  em->grid2D(grid);
  int nx = grid.size();
  int ny = grid[0].size();

  Eigen::SparseMatrix<double> K = em->getStiffnessSparse();
  fixRigid2d(K,em);
  enforcePeriodicity(K, em, grid);

  double olddepth = em->depth[0];
  em->depth[0] = 1;
  Eigen::MatrixXd Ke=em->getStiffness(0);
  em->depth[0] = olddepth;
  std::cout<<Ke<<"\n";

  linsolver.init();
  linsolver.init(K);

  em->depthColorScale = 1;
  srand(std::time(0));

  int len = microstructures[0].size();
  //index of structures with interesting property.
  std::vector<int> goodmat;
  for(unsigned int ii=0; ii<microstructures.size(); ii++){
    goodmat.push_back(ii);
  }
  //for each microstructure
  for(unsigned int gi = 0; gi<goodmat.size(); gi++){
    int mi = goodmat[gi];
    for(int jj = 0 ; jj<len; jj++){
      em->depth[jj] = microstructures[mi][jj];
    }
    clampVec(em->depth, minDensity, maxDensity);

    solveLinStatic(em, grid);
    double v0 = G(em, grid);
    double dx, dy ; 
    int idx = trajectory.size()-2;
    dx = trajectory[idx];
    dy = trajectory[idx+1]; 
    continue;
    if(dy>0){
//    if(dy/dx<-0.9){
    }else{
      trajectory.pop_back();
      trajectory.pop_back();
      continue;
    }
    dx0 = 0.25*dx;
    dy0 = 0.5*dx;
    double shrink = 0.5;
    for(unsigned int jj =0; jj<em->depth.size(); jj++){
      em->depth[jj] = 0.5 + shrink*(em->depth[jj]-0.5);
    }
//    readGrid2D(matgrid, em->depth, nx, ny);
    solveLinStatic(em, grid);
    v0 = G(em, grid);
    std::vector<double> x0;
    int logIdx = 0;
    x0 = em->depth;
    //continuous search
    for(int iter = 0; iter<ITER; iter++){
      std::vector<double> dgdx = DGDx(em, grid, Ke);
//      std::vector<double> grad = DGDp(matgrid, dgdx, nx, ny);
      std::vector<double> grad = dgdx;

      int maxidx = max_abs(grad);
      double maxval = std::abs(grad[maxidx]);
      if(maxval<1e-20){
        break;
      }
      double scale = maxStep/maxval;
      scaleVec(grad, scale);
      log<<iter<<"\n";

      //line search
      float h = 1;
      std::cout<<"v0:\n";
      double v0 = G(em, grid);
//      std::vector<double> p0 = matgrid.param;
      std::vector<double> p0 = em->depth;
      bool fail = false;
      while(1){
//        matgrid.param = p0;
//        addTimes(matgrid.param, -h, grad);
//        clampVec(matgrid.param, minDensity, maxDensity);
//        readGrid2D(matgrid, em->depth, nx,ny);

        em->depth = p0;
        addTimes(em->depth, -h, grad);
        clampVec(em->depth, minDensity, maxDensity);

        solveLinStatic(em, grid);
        double v1 = G(em, grid);
        int idx = trajectory.size()-2;
        dx = trajectory[idx];
        dy = trajectory[idx+1];
        std::cout<<"h "<<h<<" "<<v0<<" "<<v1<<"\n";
        std::cout<<"dx dy"<<dx<<" "<<dy<<"\n";
        if(v1<v0){
          break;
        }else{
          h *= 0.5;
        }
        if(h<1e-5){
          fail = true;
          break;
        }
      }
      if(fail){
        break;
      }
      for(int ii = logIdx; ii<trajectory.size(); ii++){
        log<<trajectory[ii]<<"\n";
      }
      logIdx = trajectory.size();
    }

    //  discreteSearch(em, stepper, grid);
    log.close();

    for(unsigned int ii = 0; ii<em->depth.size(); ii++){
      int val = 0;
      int xi = ii/ny;
      int yi = ii%ny;
      int idx = yi*nx+xi;
      if(em->depth[idx]>0.4){
        val = 1;
      }
    //  std::cout<<val<<" ";
      if(yi==ny-1){
    //    std::cout<<"\n";
      }
    }
  }
  out_eval<<trajectory.size()<<"\n";
  for(unsigned int ii= 0; ii<trajectory.size(); ii++){
    out_eval<<trajectory[ii]<<"\n";
  }

  out_eval.close();

  printVec(em->depth, nx, ny);
}
