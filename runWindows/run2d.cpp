#include "runtime.hpp"
#include "cfgUtilities.h"
#include "EigenUtil.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "FileUtil.hpp"
#include "PiecewiseConstant2D.hpp"
#include "PiecewiseConstantSym2D.hpp"
#include "PiecewiseConstantCubic2D.hpp"

#include "MaterialQuad2D.h"
#include "StrainLin2D.h"
#include <algorithm>
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
void loadBinary(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments);

void run2D(const ConfigFile & conf)
{
  int nx = 16;
  int ny = 16;
  Eigen::VectorXd x0;

  if (conf.hasOpt("fileRootName")){
    loadBinary(conf, continuousMaterials);
  }
  if (conf.hasOpt("structurelist")){
    loadText(conf, continuousMaterials);
  }
    
  x0 = 0.5 * Eigen::VectorXd::Ones(nx * ny);
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
  PiecewiseConstantSym2D * field = new PiecewiseConstantSym2D();
  //PiecewiseConstantCubic2D * field = new PiecewiseConstantCubic2D();
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
    return;
  }
  std::string task = conf.getString("task");
  std::thread thread;
  if (render){
    if (task == "opt"){
      thread = std::thread(optMat, fem, nSteps);
    }
    else if (task == "compute"){
      thread = std::thread(computeMat, fem, conf);
    }

    Render render;
    World * world = new World();
    world->em2d.push_back(em);
    render.init(world);
    render.loop();
  }
  else{
    if (task == "opt"){
      optMat(fem, nSteps);
    }
    else if (task == "compute"){
      computeMat(fem, conf);
    }
  }
}

void printStructure(FEM2DFun * fem, std::ostream & out)
{
  out << fem->m_nx << " " << fem->m_ny << "\n";
  for (int ii = 0; ii < fem->m_nx; ii++){
    for (int jj = 0; jj < fem->m_ny; jj++){
      double val = fem->distribution[ii*fem->m_ny + jj];
      out << val << " ";
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

void shrinkVector(Eigen::VectorXd & x, const std::vector<double> & a, float shrink_ratio)
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
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      x[2 * ii*ny + 2 * jj] = a[ii*ny + jj];
      x[2 * ii*ny + 2 * jj + 1] = a[ii*ny + jj];
      x[2 * (ii + 1)*ny + 2 * jj] = a[ii*ny + jj];
      x[2 * (ii + 1)*ny + 2 * jj + 1] = a[ii*ny + jj];
    }
  }
  return x;
}

std::string sequenceFilename(std::string prefix, int idx, std::string suffix)
{
  return prefix + std::to_string(idx) + suffix;
}

std::vector< std::vector<double> > interpArr(const std::vector<std::vector<double> > & a, int nStep, double lb, double ub)
{
  std::vector< std::vector<double> > out;
  for (unsigned int ii = 0; ii < a.size(); ii++){
    std::vector<double> b = a[ii];
    for (int jj = 0; jj < nStep; jj++){
      for (unsigned int kk = 0; kk < b.size(); kk++){
        b[kk] = (a[ii][kk] - lb) * (jj + 1.0) / nStep + lb;
      }
      out.push_back(b);
    }
    for (int jj = 0; jj < nStep; jj++){
      for (unsigned int kk = 0; kk < b.size(); kk++){
        b[kk] = a[ii][kk] + (jj + 1.0) / nStep * (ub - a[ii][kk]);
      }
      out.push_back(b);
    }
  }
  return out;
}

void fill(Eigen::VectorXd & x, double v)
{
  for (int ii = 0; ii < x.size(); ii++){
    x[ii] = v;
  }
}

void computeMat(FEM2DFun * fem, const ConfigFile & conf)
{
  RealField * field = fem->field;
  Eigen::VectorXd x1 = fem->param;
  Eigen::Vector3d fx = fem->forceMagnitude * Eigen::Vector3d(1, 0, 0);
  Eigen::Vector3d fy = fem->forceMagnitude * Eigen::Vector3d(0, 1, 0);
  int nrows = (int)fem->externalForce[0].size();
  fem->externalForce.resize(3);
  fem->externalForce[1].resize(nrows);
  fem->externalForce[2].resize(nrows);
  stretchY(fem->em, fy, fem->grid, fem->externalForce[1]);
  //shear force
  stretchY(fem->em, fx, fem->grid, fem->externalForce[2]);
  for (unsigned int ii = 0; ii < fem->em->fe.size(); ii++){
    fem->em->fe[ii][0] = fem->externalForce[2][ii*2];
    fem->em->fe[ii][1] = fem->externalForce[2][ii * 2+1];
  }
  fem->initArrays();
  std::string outputprefix("paraml");
  if (conf.hasOpt("outputprefix")){
    outputprefix = conf.getString("outputprefix");
  }
  if (conf.hasOpt("interp")){
    int nStep = conf.getInt("interp");
    for (unsigned int ii = 0; ii < continuousMaterials.size(); ii++){
      for (unsigned int jj = 0; jj < continuousMaterials[ii].size(); jj++){
        continuousMaterials[ii][jj] = std::max(1e-3, continuousMaterials[ii][jj]);
      }
    }
    continuousMaterials = interpArr(continuousMaterials, nStep, 1e-3, 1);
  }
  std::string materialOutName = sequenceFilename(outputprefix, chunkIdx, ".txt");
  materialOut.open(materialOutName);
  std::ofstream structOut("s.txt");
  for (int ii = 0; ii < continuousMaterials.size(); ii++){
    Eigen::VectorXd x0(continuousMaterials[ii].size());
    shrinkVector(x0, continuousMaterials[ii], 1.0);
    x1 = firstQuadrant(x0, fem->m_nx, fem->m_ny);
    ///check solid material
    fill(x1, 1.0);
    fem->setParam(x1);
    double dx = fem->dx;
    double dy = fem->dy;
    //stretch in y direction under stretching force in y
    double dyy = measureStretchY(fem->em, fem->u[1], fem->grid);
    double dxy = measureShearX(fem->em, fem->u[2], fem->grid);
    double F = fem->forceMagnitude;
    materialOut << F/dx << " " << F/dyy << " " << -dy/dx << " " << F/dxy<<" "<<fem->density << "\n";
    for (unsigned int jj = 0; jj < fem->distribution.size(); jj++){
      structOut << fem->distribution[jj] << " ";
    }
    structOut << "\n";
  }
  materialOut.close();
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
  std::string materialOutName = sequenceFilename("materialStructure", chunkIdx, ".txt");
  materialOut.open(materialOutName);
  double shrink_ratio = 0.3;
  for (int ii = endIdx - 2; ii < endIdx; ii++){
    Eigen::VectorXd x0(materialAssignments[ii].size());
    shrinkVector(x0, materialAssignments[ii], shrink_ratio);
    x1 = firstQuadrant(x0, fem->m_nx, fem->m_ny);
    x1 = upsampleVector(x1, 2 * fem->m_nx, 2 * fem->m_ny);
    fem->setParam(x1);
    //printStructure(fem);
    fem->dx0 = fem->dx;
    //set up objective for negative poisson ratio.
    fem->dy0 = fem->dx;
    fem->m0 = 0.3*sum(fem->distribution) / fem->distribution.size();
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

int lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h)
{
  double f0 = fun->f();
  Eigen::VectorXd grad0 = fun->df();
  double norm0 = infNorm(grad0);
  int nSteps = 10;
  int ret = -1;
  for (int step = 0; step<nSteps; step++){
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
  std::cout << "max diff " << max_diff << " " << max_val << "\n";
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
    Eigen::VectorXd fe = Eigen::Map<Eigen::VectorXd>(fem->externalForce[ii].data(), fem->externalForce[ii].size());

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

void loadText(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments)
{
  std::string filename("materialstructure.txt");
  if (conf.hasOpt("structurelist")){
    filename = conf.getString("structurelist");
  }
  FileUtilIn in(filename.c_str());
  int nMat, gridSize;
  in.in >> nMat >> gridSize;
  materialAssignments.resize(nMat);
  for (int ii = 0; ii < nMat; ii++){
    materialAssignments[ii].resize(gridSize);
    for (int jj = 0; jj < gridSize; jj++){
      in.in >> materialAssignments[ii][jj];
    }
  }
  in.close();
}

void loadBinary(const ConfigFile & conf, std::vector<std::vector<int> > & materialAssignments)
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
  std::string matAssignmentFile = conf.dir + "/" + fileRootName + fileExtension;
  success = cfgUtil::readBinary<int>(matAssignmentFile, materialAssignments);
  if (!success){
    std::cout << "Can't read " << matAssignmentFile << "\n";
  }

  std::ofstream out("text.txt");
  out << materialAssignments.size() << " " << materialAssignments[0].size() << "\n";
  for (unsigned int ii = 0; ii < materialAssignments.size(); ii++){
    for (unsigned int jj = 0; jj < materialAssignments[ii].size(); jj++){
      int val = materialAssignments[ii][jj];
      if (val == 0){
        out << 1e-3 << " ";
      }
      else{
        out << materialAssignments[ii][jj] << " ";
      }
    }
    out << "\n";
  }
  out.close();
}

void loadBinary(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments)
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
  std::string matAssignmentFile = conf.dir + "/" + fileRootName + fileExtension;
  success = cfgUtil::readBinary<double>(matAssignmentFile, materialAssignments);
  if (!success){
    std::cout << "Can't read " << matAssignmentFile << "\n";
  }
}

void transposeSave(std::vector<std::vector<double> > & continuousMaterials, int nx, int ny)
{
  //transpose each structure and save.
  unsigned int N = continuousMaterials.size();
  std::cout << N << "N\n";
  for (unsigned int si = 0; si< N; si++){
    std::vector<double> mat = continuousMaterials[si];
    std::cout << si << "\n";
    for (int ii = 0; ii < 16; ii++){
      for (int jj = 0; jj < 16; jj++){
        mat[ii * 16 + jj] = continuousMaterials[si][jj * 16 + ii];
      }
    }
    continuousMaterials.push_back(mat);
  }
  cfgUtil::writeBinary("mat.bin", continuousMaterials);
}
