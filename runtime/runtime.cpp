#include "ConfigFile.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "StepperGrad.hpp"
#include "StepperNewton.hpp"
#include "AdmmCPU.hpp"
#include "AdmmNoSpring.hpp"
//#include "IpoptStepper.hpp"
#include <thread>
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "StrainCorotLin.hpp"
#include "UnitTests.hpp"

void runSim(ElementMesh * m, Stepper * stepper)
{
  stepper->step(m);
}

void runTest()
{
  //ElementCoarseTest();
  stiffnessTest(0);
  //testCoarseDefGrad();
	//forceTest(0);
  system("pause");
  exit(0);
}

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);
  if (conf.getBool("test")){
    runTest();
  }
  int nx = 4, ny=16, nz=4;
  //int nx = 1, ny=4, nz=1;

  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  //StrainEneNeo ene;
  //std::vector<StrainCorotLin> ene(2);
  std::vector<StrainLin> ene(2);

  ene[0].param[0] = 10000;
  ene[0].param[1] = 100000;
  ene[1].param[0] = 100000;
  ene[1].param[1] = 1000000;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
    em->m.push_back(&material[ii]);
  }

  //assign some materials
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      for (int kk = 0; kk < nz; kk++){
        int eidx = em->GetEleInd(ii, jj, kk);
        if (jj < ny / 2){
        //  em->me[eidx] = 1;
        }
      }
    }
  }

  Vector3f ff(1.25,-5,0);
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
        em->fe[vidx] = ff;
      }
    }
  }
  em->check();

  World * world = new World();
  world->em.push_back(em);
  
  Stepper * stepper = 0;
  int nSteps = 10000;
  if (conf.hasOpt("nSteps")){
    nSteps = conf.getInt("nSteps");
  }
  std::string stepperType = conf.getString("stepper");
  if (stepperType == "newton"){
    stepper = new StepperNewton();
  }
  else{
    stepper = new AdmmNoSpring();
  }
  //stepper->rmRigid = true;
  //IpoptStepper * stepper = new IpoptStepper();
  //AdmmCPU * stepper = new AdmmCPU();
  //stepper->ro0 = 500;
  
  stepper->nSteps = nSteps;
  std::thread simt(runSim, em, stepper);
  Render render;
  render.init(world);
  render.loop();
  simt.join();
	return 0;
}
