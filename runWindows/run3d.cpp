#include "runtime.hpp"
#include "cfgUtilities.h"
#include "EigenUtil.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "FEM3DFun.hpp"
#include "FileUtil.hpp"
#include "PiecewiseConstant3D.hpp"
#include "PiecewiseConstantSym3D.hpp"
#include "PiecewiseConstantCubic3D.hpp"

#include "MaterialQuad.hpp"
#include "StrainLin.hpp"
#include <algorithm>
#include <iomanip>
#include <thread>

extern std::ofstream logfile;

struct Opt3DArgs
{
  FEM3DFun * fem;
  const ConfigFile * conf;
};

void shrinkVector(Eigen::VectorXd & x, float shrink_ratio)
{
  for (int ii = 0; ii < x.size(); ii++){
    x[ii] = 0.5 + shrink_ratio * (x[ii] - 0.5);
  }
}

void optMat3D(Opt3DArgs * arg)
{
  int nSteps = 500;
  const ConfigFile * conf = arg->conf;
  FEM3DFun * fem = arg->fem;

  //RealField * field = fem->field;
  
    //for (int ii = 0; ii < x1.size(); ii++){
  //  x1[ii] = 0.5 + (rand() / (float)RAND_MAX - 0.5) * 0.2;
  //}
 
  std::vector<std::vector<double> > structures;
  //std::vector<int> idx;
  //cfgUtil::readBinary<int>("../data/boundaryPoints.bin", idx);
  loadIntBinary(*conf, structures);
  std::ofstream matStruct("struct.txt");
  double shrinkRatio = 0.3;
  matStruct << fem->gridSize[0] << " " << fem->gridSize[1] << " " << fem->gridSize[2] << "\n";
  for (unsigned int si = 17; si < structures.size(); si++){
    std::cout << si << "\n";
    std::vector<double> s3d = structures[si];
    Eigen::VectorXd x1 = fem->param;
    //copy 8x8x8 corner to 3D structure.
    std::vector<int> paramSize = fem->gridSize;
    for (int dim = 0; dim < (int)paramSize.size(); dim++){
      paramSize[dim] /= 2;
    }
    for (int ix = 0; ix < paramSize[0]; ix++){
      for (int iy = 0; iy < paramSize[1]; iy++){
        for (int iz = 0; iz < paramSize[2]; iz++){         
          int linearIdx = ix * paramSize[1] * paramSize[2] + iy *paramSize[2] + iz;
          int inputIdx = 0;
          if (s3d.size() == x1.size()){
            //need to upsample s3d by a factor of 2.
            inputIdx = ix / 2 * paramSize[1] * paramSize[2] + iy / 2 * paramSize[2] + iz / 2;
          }
          else{
            inputIdx = ix * fem->gridSize[1] * fem->gridSize[2] + iy * fem->gridSize[2] + iz;
          }
          x1[2*linearIdx] = s3d[inputIdx];
          //x1[2*linearIdx+1] = 1-s3d[inputIdx];
          x1[2 * linearIdx + 1] = 1;
        }
      }
    }
    clampVector(x1, fem->lowerBounds, fem->upperBounds);
    fem->setParam(x1);
    double E = 1.0 / fem->G(0, 0);
    double nu = -fem->G(1, 0) / fem->G(0, 0);
    std::cout << fem->G(0, 0) << " " << fem->G(1, 0) << " " << fem->G(2, 0) << " " << fem->G(3, 3) << "\n";
    std::cout << "rho: " << fem->density<<", E: " << E << ", nu: " << nu << "\n";
    double val = fem->f();
    fem->m0 = 2*fem->density;
    //scale mass term to roughly displacement term.
    fem->mw = 0.1 * fem->G(0, 0) / fem->density;
    double Et = 2*E;
    nu = 2 * nu;
    fem->G0 = fem->G;
    //poisson's ratio objective
    fem->G0(0, 0) = 1.0 / Et;
    fem->G0(1, 0) = -nu * fem->G0(0, 0);// fem->G0(0, 0);
    fem->G0(2, 0) = -nu * fem->G0(0,0);
    int input;
    std::cin >> input;
    ////std::cout << fem->G << "\n";
    ////for test only look at the first displacement.
    for (int ii = 1; ii < fem->wG.size(); ii++){
      fem->wG(ii) = 0;
    }
//fem->G0(3, 3) *= 0.3;
//    fem->wG(3) = 0.2;
    shrinkVector(x1, shrinkRatio);
    ////logfile.open("log3d.txt");
    ////check_df(fem, x1, 1e-3);
    for (int j = 0; j < 50; j++){
      nSteps = 10;
      gradientDescent(fem, x1, nSteps);
      for (unsigned int ii = 0; ii < fem->distribution[0].size(); ii++){
        matStruct << fem->distribution[0][ii] << " " <<fem->distribution[1][ii]<<"\n";
        if (ii%fem->gridSize[0] == 0){
          matStruct << "\n";
        }
      }
      matStruct << "\n";
    }
  }
  matStruct.close();
  system("pause");
}

void run3D(const ConfigFile & conf)
{
  bool render = conf.getBool("render");
  std::string task = conf.getString("task");

  int gridres = 32;
  if (conf.hasOpt("gridres")){
    gridres = conf.getInt("gridres");
  }
  if (conf.hasOpt("maxStep")){
    maxStep = conf.getFloat("maxStep");
  }
  int nx = gridres;
  int ny = gridres;
  int nz = gridres;
  int nSteps =500;
  Eigen::VectorXd x0;
  
  ElementRegGrid * em = new ElementRegGrid(nx, ny, nz);
  std::vector<StrainLin> ene(1);
  //E = 1e3.
  ene[0].param[0] = 100;
  ene[0].param[1] = 2400;
  ene[0].param[1] = 900;

  std::vector<MaterialQuad * > material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    material[ii] = new MaterialQuad();
    for (unsigned int jj = 0; jj < material[ii]->e.size(); jj++){
      material[ii]->e[jj] = &ene[ii];
    }
    em->addMaterial(material[ii]);
  }
  em->initArrays();
  em->check();

  FEM3DFun * fem = new FEM3DFun();

  //PiecewiseConstant3D * field = new PiecewiseConstant3D();
  //PiecewiseConstantSym3D * field = new PiecewiseConstantSym3D();
  PiecewiseConstantCubic3D * field = new PiecewiseConstantCubic3D();
  field->allocate(nx / 2, ny / 2, nz / 2 );
  PiecewiseConstantCubic3D * field1 = new PiecewiseConstantCubic3D();
  field1->allocate(nx / 2, ny / 2, nz / 2);
  //field->allocate(nx , ny , nz);
  fem->lowerBounds = 0 * Eigen::VectorXd::Ones(2*field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(2*field->param.size());
   
  float fx = 1;
  fem->forceMagnitude = (double)fx;
  
  fem->m_fixRigid = true;
  fem->m_periodic = true;
  fem->constrained = false;

  fem->em = em;
  fem->field[0] = field;
  fem->field[1] = field1;
  fem->gridSize[0] = nx;
  fem->gridSize[1] = ny;
  fem->gridSize[2] = nz;
  fem->m0 = 0.4;
  fem->mw = 0.1;
  field->param = Eigen::VectorXd::Ones(field->param.size());
  x0 = fem->upperBounds;
  fem->init(x0);
  std::thread thread;
  if (render){
    if (task == "opt"){
      Opt3DArgs * arg = new Opt3DArgs;
      arg->fem = fem;
      arg->conf = &conf;
      thread = std::thread(optMat3D, arg);
    }
    Render render;
    World * world = new World();
    world->em.push_back(em);
    world->u = &fem->u;
    world->fe = &fem->externalForce;
    render.init(world);
    render.loop();
  }
}
