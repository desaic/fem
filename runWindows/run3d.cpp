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

  RealField * field = fem->field;
  
    //for (int ii = 0; ii < x1.size(); ii++){
  //  x1[ii] = 0.5 + (rand() / (float)RAND_MAX - 0.5) * 0.2;
  //}
 
  std::vector<std::vector<double> > structures;
  //std::vector<int> idx;
  //cfgUtil::readBinary<int>("../data/boundaryPoints.bin", idx);
  loadIntBinary(*conf, structures);
  //std::ofstream matStruct("struct.txt");
  double shrinkRatio = 0.3;
  std::vector<float> param;
  structures.push_back(std::vector<double>(4096*8, 1));
  for (unsigned int si = 0; si < structures.size(); si++){
    std::cout << si << "\n";
    std::vector<double> s3d = structures[si];
    Eigen::VectorXd x1 = fem->param;
    //copy 8x8x8 corner to 3D structure.
    for (int ii = 0; ii < x1.size(); ii++){
      x1[ii] = 1;// 1e-3;
    }
    std::vector<int> paramSize(3, 0);
    for (int dim = 0; dim < (int)paramSize.size(); dim++){
      paramSize[dim] = fem->gridSize[dim];// / 2;
    }
    for (int ix = 0; ix < paramSize[0]; ix++){
      for (int iy = 0; iy < paramSize[1]; iy++){
        for (int iz = 0; iz < paramSize[2]; iz++){
          x1[ix * paramSize[1] * paramSize[2] + iy *paramSize[2] + iz] = 0.1;// s3d[ix * fem->gridSize[1] * fem->gridSize[2] + iy * fem->gridSize[2] + iz];
          if (ix % 2 == 0){
            x1[ix * paramSize[1] * paramSize[2] + iy *paramSize[2] + iz] = 1;
          }
        }
      }
    }
    clampVector(x1, fem->lowerBounds, fem->upperBounds);
    fem->setParam(x1);
    param.push_back(fem->density);
    param.push_back(1.0 / fem->G(0, 0));
    param.push_back(-fem->G(1, 0)/fem->G(0,0));
    param.push_back(1.0 / fem->G(3, 3));
    std::cout << fem->G(0, 0) << " " << fem->G(1, 0) << " " << fem->G(2, 0) << " " << fem->G(3, 3) << "\n";
    int input;
    std::cin >> input;
    //double val = fem->f();
    //fem->m0 = 0.5 * fem->density;
    //fem->mw = 0.1 * fem->G(0, 0) / fem->density;
    //fem->G0 = fem->G;
    ////-0.45 poisson's ratio objective
    //fem->G0(1, 0) = -1 * fem->G0(0, 0);
    //fem->G0(2, 0) = -1 * fem->G0(0, 0);

    //std::cout << fem->G << "\n";
    //for test only look at the first displacement.
    //for (int ii = 1; ii < fem->wG.size(); ii++){
    //  fem->wG(ii) = 0;
    //  //  fem->wG(ii) = 1 + (rand() / (float)RAND_MAX - 0.5)*0.3 ;
    //}
    //shrinkVector(x1, shrinkRatio);
    //logfile.open("log3d.txt");
    //check_df(fem, x1, 1e-3);
    //scale mass term to roughly displacement term.
    //gradientDescent(fem, x1, nSteps);
  //  for (unsigned int ii = 0; ii < fem->distribution.size(); ii++){
  //    matStruct << fem->distribution[ii] << " ";
  //  }
  //  matStruct << "\n";
  }
  //matStruct.close();
    cfgUtil::writeBinary<float>("param2d.bin", param);
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

  PiecewiseConstant3D * field = new PiecewiseConstant3D();
  //PiecewiseConstantSym3D * field = new PiecewiseConstantSym3D();
  //PiecewiseConstantCubic3D * field = new PiecewiseConstantCubic3D();
  //field->allocate(nx / 2, ny / 2, nz / 2 );
  field->allocate(nx , ny , nz);
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());
   
  float fx = 1;
  fem->forceMagnitude = (double)fx;


  fem->m_fixRigid = false;
  fem->m_periodic = false;
  fem->constrained = true;

  fem->em = em;
  fem->field = field;
  fem->gridSize[0] = nx;
  fem->gridSize[1] = ny;
  fem->gridSize[2] = nz;
  fem->m0 = 0.4;
  fem->mw = 0.1;
  field->param = Eigen::VectorXd::Ones(field->param.size());
  x0 = field->param;
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
