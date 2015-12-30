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

#include "MaterialQuad.hpp"
#include "StrainLin.hpp"
#include <algorithm>
#include <iomanip>
#include <thread>

void computeMat3D(FEM3DFun * fem, const ConfigFile & conf)
{

}

void optMat3D(FEM3DFun * fem, int nSteps)
{
  RealField * field = fem->field;
  Eigen::VectorXd x1 = fem->param;
  fem->setParam(x1);
  fem->m0 = sum(fem->distribution) / fem->distribution.size();
  //scale mass term to roughly displacement term.
  gradientDescent(fem, x1, nSteps);
}

void run3D(const ConfigFile & conf)
{
  bool render = conf.getBool("render");
  std::string task = conf.getString("task");
  int nx = 16;
  int ny = 16;
  int nz = 16;
  int nSteps = 10;
  Eigen::VectorXd x0;
  x0 =  Eigen::VectorXd::Ones(nx * ny * nz);
  ElementRegGrid * em = new ElementRegGrid(nx, ny, nz);
  std::vector<StrainLin> ene(1);
  ene[0].param[0] = 3448.275862;
  ene[0].param[1] = 31034.48276;

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
  field->allocate(nx , ny , nz );

  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());

  fem->em = em;
  fem->field = field;
  fem->gridSize[0] = nx;
  fem->gridSize[1] = ny;
  fem->gridSize[2] = nz;
  fem->m0 = 0.4;
  fem->mw = 0.1;

  fem->init(x0);
  std::thread thread;
  if (render){
    if (task == "opt"){
      thread = std::thread(optMat3D, fem, nSteps);
    }
    else if (task == "compute"){
      thread = std::thread(computeMat3D, fem, conf);
    }

    Render render;
    World * world = new World();
    world->em.push_back(em);
    render.init(world);
    render.loop();
  }
}
