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
#include "UnitTests.hpp"

void runSim(ElementMesh * m, Stepper * stepper)
{
  stepper->step(m);
}

void runTest()
{
  //ElementCoarseTest();
  stiffnessTest();
	//forceTest();
  system("pause");
  exit(0);
}

int main(int argc, char* argv[])
{
 // runTest();
  int nx = 2, ny=8, nz=2;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  //StrainEneNeo ene;
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
          em->me[eidx] = 1;
        }
      }
    }
  }

  Vector3f ff(5,-20,0);
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
  
  //StepperNewton *stepper= new StepperNewton();
  //stepper->rmRigid = true;
  //IpoptStepper * stepper = new IpoptStepper();
  AdmmCPU * stepper = new AdmmCPU();
  stepper->ro0 = 500;
  //AdmmNoSpring *stepper = new AdmmNoSpring();
  stepper->nSteps = 1000;
  std::thread simt(runSim, em, stepper);
  Render render;
  render.init(world);
  render.loop();
  simt.join();
	return 0;
}
