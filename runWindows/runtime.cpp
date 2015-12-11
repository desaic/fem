#include "ConfigFile.hpp"
#include "EigenUtil.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "StepperGrad.hpp"
#include "StepperNewton.hpp"
//#include "IpoptStepper.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "StrainCorotLin.hpp"
#include "UnitTests.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

//reorder stiffness matrix for a single element to compare to vega
void OrderVega(Eigen::MatrixXd & K){
  if (K.rows() != 24){
    //only works with 24x24 matrix
    return;
  }
  Eigen::MatrixXd Kc = K;
  int vidx[8] = { 0, 4, 6, 2, 1, 5, 7, 3 };
  for (int ii = 0; ii < 8; ii++){
    for (int jj = 0; jj < 8; jj++){
      K.block(3 * ii, 3 * jj, 3, 3) = Kc.block(3 * vidx[ii], 3 * vidx[jj], 3, 3);
    }
  }
}

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);  
  
  int refine = 1;
  if (conf.hasOpt("refine")){
    refine = conf.getInt("refine");
  }
  int res = (int)std::pow(2, refine);
  int nx = res, ny=4*res, nz=res;
  //int nx = 32, ny = 80, nz = 32;
  //int nx = 16, ny = 40, nz = 16;
  Vector3f ff(5, -100, 0);
  //per element pushing force
  ff = (1.0f / (nx*nz)) * ff;

  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  std::vector<StrainEneNeo> ene(2);
  //std::vector<StrainCorotLin> ene(2);
  //std::vector<StrainLin> ene(2);

  ene[0].param[0] = 34482.75862;
  ene[0].param[1] = 310344.8276;
  ene[1].param[0] = 100000;
  ene[1].param[1] = 1000000;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
    em->addMaterial(&material[ii]);
  }
  em->initArrays();
  MatrixXf K = em->getStiffness(0);
  Eigen::MatrixXd Ke(K.mm, K.nn);
  for (int ii = 0; ii < K.mm; ii++){
    for (int jj = 0; jj < K.nn; jj++){
      Ke(ii, jj) = K(ii, jj);
    }
  }
  OrderVega(Ke);
  
  std::ofstream out("K.txt");
  write_vega_lists(out, Ke.sparseView());
  out.close();

  for(int ii = 0;ii<nx;ii++){
    for(int jj =0;jj<nz;jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      for(int kk = 0;kk<4;kk++){
        int vidx =em->e[eidx]->at(aa[kk]);
        em->fixed[vidx] = 1;
      }

      eidx= em->GetEleInd(ii,ny-1,jj);
      int bb[4] = {2,3,6,7};
      for(int kk = 0;kk<4;kk++){
        int vidx = em->e[eidx]->at(bb[kk]);
        em->fe[vidx] += ff;
      }
    }
  }
  em->check();

  World * world = new World();
  bool renderOpt = true;
  if (conf.hasOpt("render") && conf.getBool("render") == false){
    renderOpt = false;
  }
  if (renderOpt){
    world->em.push_back(em);
  }
  Stepper * stepper = 0;
  int nSteps = 10000;
  if (conf.hasOpt("nSteps")){
    nSteps = conf.getInt("nSteps");
  }
  
  stepper = new StepperNewton();

  //stepper->rmRigid = true;
  
  stepper->nSteps = nSteps;
  stepper->init(em);
  stepper->launchThread();
  Render render;
  world->stepper = stepper;
  render.init(world);
  render.loop();
	return 0;
}
