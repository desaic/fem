#include "ConfigFile.hpp"
#include "EigenUtil.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "FEM2DFun.hpp"
#include "PiecewiseConstant2D.hpp"

#include "MaterialQuad2D.h"
#include "StrainLin2D.h"

#include <fstream>
#include <iostream>
#include <iomanip>

///@brief verify consistency of f and df using central differencing
void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h);

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
  PiecewiseConstant2D * field = new PiecewiseConstant2D();
  fem->em = em;
  fem->field = field;
  fem->m_nx = nx;
  fem->m_ny = ny;
  Eigen::VectorXd x0 = 0.5 * Eigen::VectorXd::Ones(em->e.size());
  fem->init(x0);
  for (int ii = 0; ii < x0.rows(); ii++){
    x0[ii] += 0.1 * (rand() / (float)RAND_MAX - 0.5);
  }
  double h = 1e-2;
  check_df(fem, x0, h);
  system("PAUSE");
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