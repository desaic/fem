#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "StepperGrad.hpp"
#include "StepperNewton.hpp"
#include "AdmmCPU.hpp"
//#include "IpoptStepper.hpp"
#include <thread>
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"

#include "UnitTests.hpp"

void runSim(ElementMesh * m, Stepper * stepper)
{
  stepper->step(m);
}

void runTest()
{
  ElementCoarseTest();
//  system("pause");
  exit(0);
}

int main(int argc, char* argv[])
{
//  runTest();
  int nx = 4, ny=16, nz=4;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  StrainEneNeo ene;
  ene.param[0] = 34482;
  ene.param[1] = 310334;
  MaterialQuad material(&ene);
  em->m.push_back(&material);
  Vector3f ff(0,-1000,0);
  for(int ii = 0;ii<nx/2;ii++){
    for(int jj =0;jj<nz/2;jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      for(int kk = 0;kk<4;kk++){
        int vidx =em->e[eidx]->at(aa[kk]);
//        em->fixed[vidx] = 1;
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
  
  StepperNewton *stepper= new StepperNewton();
  stepper->rmRigid = true;
//  IpoptStepper * stepper = new IpoptStepper();
  //AdmmCPU *stepper= new AdmmCPU();
  //stepper->ro0 = 1000;
  stepper->nSteps = 100000;
  std::thread simt(runSim, em, stepper);
  Render render;
  render.init(world);
  render.loop();
  simt.join();
	return 0;
}
