#include "ConfigFile.hpp"
#include "EigenUtil.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "FEM2DFun.hpp"
#include "PiecewiseConstant2D.hpp"
#include "PiecewiseConstantSym2D.hpp"

#include "MaterialQuad2D.h"
#include "StrainLin2D.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

///@brief sparse matrix vector product.
Eigen::VectorXd
spmv(const std::vector<int> & I, const std::vector<int> & J, const std::vector<double> & val,
const std::vector<double> x, bool triangle= true);

///@brief verify consistency of f and df using central differencing
void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h);

///@brief verify that linear static simulation is working.
///Checks K*u - f.
///Also checks df/dx = K.
///Numerical differencing is very slow. It computes all N*N entries in K, where N = Number of elements.
///@TODO: Move code to UnitTests.
void check_sim(FEM2DFun * fem, const Eigen::VectorXd & x);

///@brief finds a step size the decreases value of fun.
///@return <0 if cannot find a positive step size.
///@param x0 starting point.
///@param dir search direction.
///@param h appropriate step size if there is one. Always positive.
int lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h);

///@param nSteps maximum number of steps.
void gradientDescent(RealFun * fun, Eigen::VectorXd & x0, int nSteps);

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);  
  
  int nx = 16, ny = 16;
  Vector3f ff(100, 0, 0);
  
  ElementRegGrid2D * em = new ElementRegGrid2D(nx,ny);
  std::vector<StrainLin2D> ene(1);

  ene[0].param[0] = 34482.75862;
  ene[0].param[1] = 310344.8276;

  std::vector<MaterialQuad2D * > material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    material[ii] = new MaterialQuad2D();
    for (unsigned int jj = 0; jj < material[ii]->e.size(); jj++){
      material[ii]->e[jj] = &ene[ii];
    }
    em->addMaterial(material[ii]);
  }
  em->initArrays();
  em->check();

  FEM2DFun * fem = new FEM2DFun();
  PiecewiseConstantSym2D * field = new PiecewiseConstantSym2D();
  //for example, the parameterization can have lower resolution. Effectively refining each cell to 2x2 block.
  field->allocate(nx/2, ny/2);
  fem->em = em;
  fem->dx0 = 2e-3;
  fem->dy0 = 2e-3;
  fem->field = field;
  fem->m_nx = nx;
  fem->m_ny = ny;
  Eigen::VectorXd x0 = 0.5 * Eigen::VectorXd::Ones(field->param.size());
  //uniform lower and upper bounds for variables. 
  //Can change for more complex material distribution scheme.
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(x0.size());
  fem->upperBounds = Eigen::VectorXd::Ones(x0.size());

  fem->init(x0);
  for (int ii = 0; ii < x0.rows(); ii++){
    x0[ii] += 0.1 * (rand() / (float)RAND_MAX - 0.5);
  }
  double h = 1e-2;
  //check_df(fem, x0, h);
  check_sim(fem, x0);
  fem->setParam(x0);
  Eigen::VectorXd grad = fem->df();
  h = 1;
  gradientDescent(fem, x0, 100);
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      std::cout << fem->distribution[ii*ny + jj]<<" ";
    }
    std::cout << "\n";
  }
  system("PAUSE");
  return 0;
}

double infNorm(const Eigen::VectorXd & a)
{
  double maxval = 0;
  for (int ii = 0; ii < a.rows(); ii++){
    maxval = std::max(maxval, std::abs(a[ii]));
  }
  return maxval;
}

void gradientDescent(RealFun * fun, Eigen::VectorXd & x0, int nSteps)
{
  //maximum movement in any parameter.
  double maxStep = 0.1;
  Eigen::VectorXd x = x0;
  for (int ii = 0; ii < nSteps; ii++){
    fun->setParam(x);
    Eigen::VectorXd grad = fun->df();
    double h = 1;
    double norm = infNorm(grad);
    h = maxStep / norm;
    int ret = lineSearch(fun, x, grad, h);
    //std::cout << ii << " " << fun->f() << " " << h << "\n";
    if (ret < 0){
      break;
    }
    x0 = x;
  }
}

int lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h )
{
  double f0 = fun->f();
  Eigen::VectorXd grad0 = fun->df();
  double norm0 = infNorm(grad0);
  while (1){
    //minus sign here if we want to minimize function value.
    Eigen::VectorXd x = x0 - h*dir;
    clampVector(x, fun->lowerBounds, fun->upperBounds);
    fun->setParam(x);
    double f1 = fun->f();
    Eigen::VectorXd grad1 = fun->df();
    double norm1 = infNorm(grad1);
    //if gradient norm or function value decrease, return.
    if (norm1 < norm0 || f1 < f0){
      x0 = x;
      break;
    }
    else{
      h /= 2;
    }
  }
  return 0;
}

void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h)
{
  fun->setParam(x0);
  Eigen::VectorXd x = x0;
  Eigen::VectorXd ana_df = fun->df();
  Eigen::VectorXd num_df = Eigen::VectorXd::Zero(x0.size());

  double max_diff = 0;
  for (int ii = 0; ii < x0.rows(); ii++){
    x[ii] += h;
    fun->setParam(x);
    double f_plus = fun->f();

    x[ii] -= 2 * h;
    fun->setParam(x);
    double f_minus = fun->f();

    x[ii] += h;

    num_df[ii] = (f_plus - f_minus) / (2 * h);
    std::cout << ii << " " << ana_df[ii] << " " << num_df[ii] << "\n";
    double diff = ana_df[ii] - num_df[ii];
    max_diff = std::max(max_diff, std::abs(diff));
  }
  std::cout << "max diff " << max_diff << "\n";
}

///@brief sparse matrix vector mult
Eigen::VectorXd
spmv(const std::vector<int> & I, const std::vector<int> & J, const std::vector<double> & val,
const std::vector<double> x, bool triangle)
{
  int nrows = (int)I.size() - 1;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(nrows);
  for (int ii = 0; ii < nrows; ii++){
    for (int jj = I[ii]; jj < I[ii + 1]; jj++){
      int idx = J[jj];
      double value = val[jj];
      b[ii] += value * x[idx];
    }
  }
  if (triangle){
    for (int ii = 0; ii < nrows; ii++){
      for (int jj = I[ii]; jj < I[ii + 1]; jj++){
        int idx = J[jj];
        if (ii == idx){
          continue;
        }
        double value = val[jj];
        b[idx] += value * x[ii];
      }
    }
  }
  return b;
}

void check_sim(FEM2DFun * fem, const Eigen::VectorXd & x)
{
  fem->setParam(x);
  for (unsigned int ii = 0; ii < fem->u.size(); ii++){
    Eigen::VectorXd prod = spmv(fem->m_I, fem->m_J, fem->m_val, fem->u[ii]);
    Eigen::VectorXd fe = Eigen::Map<Eigen::VectorXd> (fem->externalForce[ii].data(), fem->externalForce[ii].size());
    
    for (int ii = 0; ii < prod.size(); ii++){
      std::cout << prod[ii] << " ";
    }
    prod -= fe; 
    double maxErr = infNorm(prod);
    std::cout << "Lin solve residual: " << maxErr << "\n";
  }
  
}
