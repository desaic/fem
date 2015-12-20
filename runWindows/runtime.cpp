#include "cfgUtilities.h"
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
#include "PiecewiseConstantCubic2D.hpp"

#include "MaterialQuad2D.h"
#include "StrainLin2D.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <thread>

//maximum movement in any parameter.
double maxStep = 0.5;
int logInterval = 10;
std::ofstream logfile;
std::ofstream materialOut;
int nChunk = 6;
int chunkIdx = 0;

std::vector<std::vector<int> > materialAssignments;
std::vector<std::vector<double> > continuousMaterials;

///@brief sparse matrix vector product.
Eigen::VectorXd
spmv(const std::vector<int> & I, const std::vector<int> & J, const std::vector<double> & val,
const std::vector<double> x, bool triangle= true);

///@brief verify consistency of f and df using central differencing
void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h);

void loadArr2D(std::vector<std::vector< int> > &arr, std::string filename);

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

void loadBinary(ConfigFile & conf, std::vector<std::vector<int> > & materialAssignments);

void printStructure(FEM2DFun * fem, std::ostream & out)
{
  out << fem->m_nx << " " << fem->m_ny << "\n";
  for (int ii = 0; ii < fem->m_nx; ii++){
    for (int jj = 0; jj < fem->m_ny; jj++){
      double val = fem->distribution[ii*fem->m_ny + jj];
      out <<  val << " ";
    }
    std::cout << "\n";
  }
}

void shrinkVector(Eigen::VectorXd & x, const std::vector<int> & a, float shrink_ratio)
{
  for (int ii = 0; ii < a.size(); ii++){
    double val = a[ii];
    val = 0.5 + shrink_ratio * (val - 0.5);
    x[ii] = val;
  }
}

Eigen::VectorXd firstQuadrant(const Eigen::VectorXd a, int nx, int ny)
{
  Eigen::VectorXd x(nx*ny / 4);
  for (int ii = 0; ii < nx / 2; ii++){
    for (int jj = 0; jj < ny / 2; jj++){
      x[ii * (ny / 2) + jj] = a[ii*ny + jj];
    }
  }
  return x;
}

Eigen::VectorXd upsampleVector(const Eigen::VectorXd a, int nx, int ny)
{
  Eigen::VectorXd x(nx * ny * 4);
  for (int ii = 0; ii < nx ; ii++){
    for (int jj = 0; jj < ny ; jj++){
      x[2 * ii*ny + 2 * jj]       = a[ii*ny + jj];
      x[2 * ii*ny + 2 * jj + 1]   = a[ii*ny + jj];
      x[2 * (ii + 1)*ny + 2 * jj] = a[ii*ny + jj];
      x[2 * (ii + 1)*ny + 2 * jj + 1] = a[ii*ny + jj];
    }
  }
  return x;
}

double sum(const Eigen::VectorXd & x)
{
  double sum = 0;
  for (int ii = 0; ii < x.size(); ii++){
    sum += x[ii];
  }
  return sum;
}

std::string sequenceFilename(std::string prefix, int idx, std::string suffix)
{
  return prefix + std::to_string(idx) + suffix;
}

//uses orthotropic structure
void optMat(FEM2DFun * fem, int nSteps)
{
  RealField * field = fem->field;
  Eigen::VectorXd x1 = fem->param;
  
  int nStructure = materialAssignments.size();
  nStructure /= nChunk;
  int startIdx = chunkIdx * nStructure;
  int endIdx = (chunkIdx + 1) * nStructure;
  //the last batch should compute untill the end.
  if (chunkIdx == nChunk){
    endIdx = std::max(endIdx, (int)materialAssignments.size());
  }
  std::string materialOutName =  sequenceFilename("materialStructure", chunkIdx, ".txt");
  materialOut.open(materialOutName);
  double shrink_ratio = 0.3;
  for (int ii = endIdx-2; ii < endIdx; ii++){
    Eigen::VectorXd x0(materialAssignments[ii].size());
    shrinkVector(x0, materialAssignments[ii], shrink_ratio);
    x1 = firstQuadrant(x0, fem->m_nx, fem->m_ny);
    x1 = upsampleVector(x1, 2 * fem->m_nx, 2*fem->m_ny);
    fem->setParam(x1);
    //printStructure(fem);
    fem->dx0 = fem->dx;
    //set up objective for negative poisson ratio.
    fem->dy0 = fem->dx;
    fem->m0 =  0.3*sum(fem->distribution)/fem->distribution.size();
    //scale mass term to roughly displacement term.
    fem->mw = 0.1 * (fem->dx0 / fem->m0);
    std::string filename;
    filename = sequenceFilename("log", ii, ".txt");
    logfile.open(filename);
    gradientDescent(fem, x1, nSteps);
    for (unsigned int jj = 0; jj < fem->distribution.size(); jj++){
      materialOut << fem->distribution[jj] << " ";
    }
    materialOut << "\n";
    logfile.close();
  }
  materialOut.close();
}

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);  
  
  int nx = 32;
  int ny = 32;
  Eigen::VectorXd x0;

  loadBinary(conf, materialAssignments);
  if (conf.hasOpt("structure")){
    std::vector<int> arr;
    std::vector<std::vector<int> > structure;
    loadArr2D(structure, conf.getString("structure"));
    nx = structure.size();
    ny = structure[0].size();
    x0 = Eigen::VectorXd::Zero(nx * ny);
    arr.resize(nx * ny);
    for (int ii = 0; ii < nx; ii++){
      for (int jj = 0; jj < ny; jj++){
        double val = structure[ii][jj];
        double shrink = 0.3;
        val = 0.5 + shrink * (val - 0.5);
        x0[ii * ny + jj] = val;
        arr[ii*ny + jj] = structure[ii][jj];
      }
    }
    materialAssignments.push_back(arr);
  } else {
    x0 = 0.5 * Eigen::VectorXd::Ones(nx * ny);
  }

  ElementRegGrid2D * em = new ElementRegGrid2D(nx, ny);
  std::vector<StrainLin2D> ene(1);
  ene[0].param[0] = 3448.275862;
  ene[0].param[1] = 31034.48276;

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
  //PiecewiseConstantSym2D * field = new PiecewiseConstantSym2D();
  PiecewiseConstantCubic2D * field = new PiecewiseConstantCubic2D();
  field->allocate(nx / 2, ny / 2);
  //for example, the parameterization can have lower resolution. Effectively refining each cell to 2x2 block.
  //PiecewiseConstant2D * field = new PiecewiseConstant2D();
  //field->allocate(nx, ny);

  //uniform lower and upper bounds for variables. 
  //Can change for more complex material distribution scheme.
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());

  fem->em = em;
  fem->dx0 = 0.22;
  fem->dy0 = 0.22;
  fem->dxw = 2;
  fem->field = field;
  fem->m_nx = nx;
  fem->m_ny = ny;
  fem->m0 = 0.4;
  fem->mw = 0.1;
  
  fem->init(x0);

  chunkIdx = conf.getInt("chunkIdx");
  int nSteps = conf.getInt("nSteps");
  bool render = conf.getBool("render");
  bool test = conf.getBool("test");
  if (conf.hasOpt("maxStep")){
    maxStep = conf.getFloat("maxStep");
  }
  if (test){
    double h = 0.001;
    nSteps = 1;
    check_df(fem, x0, h);
    //check_sim(fem, x0);
    return 0;
  }

  if (render){
    std::thread thread(optMat, fem, nSteps);
    Render render;
    World * world = new World();
    world->em2d.push_back(em);
    render.init(world);
    render.loop();
  }
  else{
    optMat(fem, nSteps);
  }
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

  Eigen::VectorXd x = x0;
  for (int ii = 0; ii < nSteps; ii++){
    fun->setParam(x);
    Eigen::VectorXd grad = fun->df();
    double h = 1;
    double norm = infNorm(grad);
    h = maxStep / norm;
    int ret = lineSearch(fun, x, grad, h);
    std::cout << ii << " " << fun->f() << " " << h << "\n";
    if (ret < 0){
      break;
    }
    x0 = x;
    if (ii % logInterval == 0){
      fun->log(logfile);
    }
  }
  fun->log(logfile);
}

int lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h )
{
  double f0 = fun->f();
  Eigen::VectorXd grad0 = fun->df();
  double norm0 = infNorm(grad0);
  int nSteps = 10;
  int ret = -1;
  for(int step = 0; step<nSteps; step++){
    //minus sign here if we want to minimize function value.
    Eigen::VectorXd x = x0 - h*dir;
    clampVector(x, fun->lowerBounds, fun->upperBounds);
    fun->setParam(x);
    double f1 = fun->f();
    Eigen::VectorXd grad1 = fun->df();
    double norm1 = infNorm(grad1);
    //if gradient norm or function value decrease, return.
    std::cout << h << " " << norm1 << " " << f1 << "\n";
    if (norm1 < norm0 || f1 < f0){
      x0 = x;
      ret = 0;
      break;
    }
    else{
      h /= 2;
    }
  }
  return ret;
}

void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h)
{
  fun->setParam(x0);
  Eigen::VectorXd x = x0;
  Eigen::VectorXd ana_df = fun->df();
  Eigen::VectorXd num_df = Eigen::VectorXd::Zero(x0.size());

  double max_diff = 0;
  double max_val = 0;
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
    max_val = std::max(max_val, std::abs(ana_df[ii]));
  }
  std::cout << "max diff " << max_diff << " "<<max_val<<"\n";
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

void loadArr2D(std::vector<std::vector< int> > &arr, std::string filename)
{
  std::ifstream in(filename);
  int nx = 0, ny = 0;
  in >> nx;
  in >> ny;
  arr.resize(nx);
  for (int ii = 0; ii < nx; ii++){
    arr[ii].resize(ny);
  }
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      in >> arr[jj][ii];
    }
  }
  in.close();
}

void loadBinary(ConfigFile & conf, std::vector<std::vector<int> > & materialAssignments)
{
  std::string fileRootName("../data/level16");
  std::string fileExtension(".bin");
  if (conf.hasOpt("fileRootName")){
    fileRootName = conf.getString("fileRootName");
  }
  if (conf.hasOpt("fileExtension")){
    fileExtension = conf.getString("fileExtension");
  }
  bool success;
  std::string matAssignmentFile = conf.dir + "/" + fileRootName + "_matAssignments" + fileExtension;
  success = cfgUtil::readBinary<int>(matAssignmentFile, materialAssignments);
  if (!success){
    std::cout << "Can't read " << matAssignmentFile << "\n";
  }

}
