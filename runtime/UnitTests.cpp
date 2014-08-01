#include "UnitTests.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include <iostream>
void forceTest()
{
  int nx = 1,ny=1,nz=1;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  StrainEneNeo ene;
  ene.param[0] = 34400;
  ene.param[1] = 310000;
  MaterialQuad material(&ene);
  em->m.push_back(&material);
  em->me.resize(em->e.size(),0);
  em->fe.resize(em->x.size());
  em->fixed.resize(em->x.size());
  Vector3f ff(0,10000,0);
  for(int ii = 0;ii<nx;ii++){
    for(int jj =0;jj<nz;jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      int bb[4] = {2,3,6,7};
      for(int kk = 0;kk<4;kk++){
        int vidx =em->e[eidx]->at(aa[kk]);
  //      em->fixed[vidx] = 1;
        vidx = em->e[eidx]->at(bb[kk]);
//        em->fe[vidx] = ff;
      }
    }
  }
  em->check();

  em->x[0][0] += 0.2f;
  em->x[1][2] -= 0.3f;
  float h = 0.0001f;
    std::vector<Vector3f>force = em->getForce();
  for(size_t ii = 0;ii<em->x.size();ii++){
    for(int jj = 0; jj<3; jj++){
      em->x[ii][jj] -= h;
      double Eminus = em->getEnergy();
      em->x[ii][jj] += 2*h;
      double Eplus = em->getEnergy();
      std::cout<<"Energy diff: "<<Eplus-Eminus<<"\n";
      double numF = (Eplus - Eminus )/(2*h);
      std::cout<<"Analytic derivative:\n"<<-force[ii][jj];
      std::cout<<"\nCentral diff:\n"<<numF<<"\n--------\n";
    }
  }
}