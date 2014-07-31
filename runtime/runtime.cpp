#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "StepperGrad.hpp"
#include <thread>
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"

void runSim(ElementMesh * m, Stepper * stepper)
{
  stepper->step(m);
}

int main(int argc, char* argv[])
{
  int nx = 2,ny=3,nz=4;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  StrainEneNeo ene;
  ene.param[0] = 34400;
  ene.param[1] = 310000;
  MaterialQuad material(&ene);
  em->m.push_back(&material);
  em->me.resize(em->e.size(),0);
  em->fe.resize(em->x.size());
  em->fixed.resize(em->x.size());
  for(int ii = 0;ii<nx;ii++){
    for(int jj =0;jj<nz;jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      for(int kk = 0;kk<4;kk++){
        em->fixed[em->e[eidx]->at(aa[kk])] = 1;
      }
    }
  }
  em->check();

  World * world = new World();
  world->em.push_back(em);
  
  StepperGrad stepper;

  std::thread simt(runSim, em, &stepper);

//  Render render;
  //render.init(world);
 // render.loop();
	return 0;
}
