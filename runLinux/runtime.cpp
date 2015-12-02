#include "ConfigFile.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementHier.hpp"
#include "ElementRegGrid.hpp"
#include "ElementMeshHier.hpp"
#include "StepperGrad.hpp"
#include "StepperNewton.hpp"
#include "StepperNewtonDyn.hpp"
//#include "NewtonCuda.hpp"
//#include "LinSolveCusp.hpp"
//#include "NewtonCuda.hpp"
#include "IpoptStepper.hpp"
#include "IpoptDynStepper.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "StrainCorotLin.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  if(argc>1){
      filename = argv[1];
  }
  ConfigFile conf;
  conf.load(filename);
  
  //make materials
  std::vector<StrainEneNeo> ene(2);
  //std::vector<StrainCorotLin> ene(2);
  //std::vector<StrainLin> ene(2);
  ene[0].param[0] = 1e6;
  ene[0].param[1] = 1e7;
  ene[1].param[0] = 1e9;
  ene[1].param[1] = 1e10;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
  }

  ElementMesh * em = 0;
  //inch to meter
  float meshScale = 0.0254;
  if(conf.hasOpt("meshfile")){
    std::string meshfile = conf.getString("meshfile");
    std::ifstream in(meshfile);
    if(!in.good()){
      std::cout<<"Can't open fem mesh "<<meshfile<<"\n";
      return -1;
    }
    em = new ElementMesh();
    em->load(in, meshScale);

    //apply some forces
    float max = 0;
    for(unsigned int ii =0 ; ii<em->x.size(); ii++){
      if(em->X[ii][1]>max){
        max = em->X[ii][1];
      }
    }
    float eps = 1e-4;
    for(unsigned int ii =0 ; ii<em->x.size(); ii++){
      if(em->X[ii][1]>max-eps && em->X[ii][0]<0.5){
        em->fe[ii] = Vector3f(0,-0.007,0);
      }
    }

    in.close();
  }else{
    int refine = 1;
    if (conf.hasOpt("refine")){
      refine = conf.getInt("refine");
    }
    int res = (int)std::pow(2, refine);
    int nx = res, ny=4*res, nz=res;
    //int nx = 32, ny = 80, nz = 32;
    //int nx = 16, ny = 40, nz = 16;
    Vector3f ff(1000, 0, 0);
    //per element pushing force
    ff = (1.0f / (nx*nz)) * ff;
    ElementRegGrid * grid = new ElementRegGrid(nx,ny,nz);
    em = grid;

    //assign some materials
    for (int ii = 0; ii < nx; ii++){
      for (int jj = 0; jj < ny; jj++){
        for (int kk = 0; kk < nz; kk++){
          int eidx = grid->GetEleInd(ii, jj, kk);
          //if (jj < ny / 2){
          //  em->me[eidx] = 1;
          //}
          //if ( (ii == 0) || ii == nx - 1 || kk == 0 || kk == nz - 1){
          //  em->me[eidx] = 1;
          //}
        }
      }
    }

    //add force and constraints
    for(int ii = 0;ii<nx;ii++){
      for(int jj =0;jj<nz;jj++){
        int eidx= grid->GetEleInd(ii,0,jj);
        int aa[4] = {0,1,4,5};
        for(int kk = 0;kk<4;kk++){
          int vidx =em->e[eidx]->at(aa[kk]);
//          em->fixed[vidx] = 1;
        }

        eidx= grid->GetEleInd(ii,ny-1,jj);
        int bb[4] = {2,3,6,7};
        for(int kk = 0;kk<4;kk++){
          int vidx = em->e[eidx]->at(bb[kk]);
          em->fe[vidx] += ff;
        }
      }
    }
  }

  for(int ii =0 ; ii<material.size(); ii++){
    em->addMaterial(&material[ii]);
  }

  em->dt = 0.01;
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
  std::string stepperType = conf.getString("stepper");
  if (stepperType == "newton"){
    stepper = new StepperNewton();
  }else if (stepperType == "newtonDyn"){
    stepper = new StepperNewtonDyn();
  }
  else if (stepperType == "newtonCuda"){
//    stepper = new NewtonCuda();
  }
  else if (stepperType == "ipopt"){
    stepper = new IpoptStepper();
  }
  else if (stepperType == "ipoptDyn"){
    stepper = new IpoptDynStepper();
  }
  else if (stepperType == "grad"){
    stepper = new StepperGrad();
  }
  else{
    //stepper = new AdmmNoSpring();
//    stepper = new ADMMStepper();
  }
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
