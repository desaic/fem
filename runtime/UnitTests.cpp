#include "UnitTests.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "MatrixXd.hpp"
#include <iostream>

void stiffnessTest()
{
  float h = 0.001f;
  ElementRegGrid * grid = new ElementRegGrid(1, 1, 1);

  StrainEneNeo* ene=new StrainEneNeo();
  ene->param[1]= 1000;
  ene->param[0] = 10000;
  MaterialQuad * material = new MaterialQuad(ene);
  grid->m.push_back(material);
  //  grid->x[1][2] +=0.01;
  //  grid->x[3][1] +=0.01;
  MatrixXd KAna = grid->getStiffness();
  int nVar = (int)grid->X.size();
  MatrixXd K(3*nVar,3*nVar);
  //check each partial derivative
  for(unsigned int ii = 0;ii<grid->x.size();ii++){
    for(int jj = 0; jj<3; jj++){
      grid->x[ii][jj] -= h;
      std::vector<Vector3f>fminus = grid->getForce();
      grid->x[ii][jj] += 2*h;
      std::vector<Vector3f>fplus = grid->getForce();
      grid->x[ii][jj] -=h;
      for(unsigned int kk = 0;kk<fminus.size();kk++){
        fplus[kk] -= fminus[kk];
        fplus[kk] /= 2*h;
        for(int ll = 0;ll<3;ll++){
          K(3*kk+ll,3*ii+jj)=-fplus[kk][ll];
        }
      }
    }
  }
  std::cout<<"Ana K:\n";
  KAna.print(std::cout);
  std::cout<<"\n";
  std::cout<<"Num K:\n";
  K.print(std::cout);
  std::cout<<"\n";

  float maxErr = 0;
  for(int ii = 0;ii<K.mm;ii++){
    for(int jj =0 ;jj<K.nn;jj++){
      float err = (float)std::abs(KAna(ii,jj)-K(ii,jj));
      if(err>maxErr){
        maxErr = err;
      }
    }
  }

  std::cout<<"max err "<<maxErr<<"\n";
}

void forceTest()
{
  int nx = 1,ny=1,nz=1;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  StrainEneNeo ene;
  ene.param[0] = 1000;
  ene.param[1] = 10000;
  MaterialQuad material(&ene);
  em->m.push_back(&material);
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
  // em->x[1][2] -= 0.3f;
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