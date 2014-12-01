#include "UnitTests.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementCoarse.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "MatrixXd.hpp"
//#include "ConjugateGradientCuda.hpp"
#include <iostream>

void ElementCoarseTest()
{
  ElementRegGrid grid (2,2,2);
  std::vector<ElementCoarse*> ce;
  ce = coarsen(grid);
  for(int ii = 0;ii<8;ii++){
    std::cout<<ce[0]->fineEle[ii];
  }

}

void cudaLinTest()
{
  //make a poisson problem
  int N = 10000;
  int nz = 3*N-4;
  std::vector<int> I(N+1),J(nz);
  std::vector<float> val, x(N,0), rhs(N,1);
  I[0] = 0;
  J[0] = 0;
  val.push_back(1);
  int nEntry=1;
  for(int ii = 1;ii<N-1;ii++){
    
    I[ii] = nEntry;
    J[nEntry]=ii-1;
    J[nEntry+1]=ii;
    J[nEntry+2]=ii+1;

    val.push_back(-1);
    val.push_back(2.15f);
    val.push_back(-1);

    nEntry+=3;
  }
  I[N-1]=nEntry;
  J[nEntry] = N-1;
  val.push_back(1);
  val[1] = 0;
  val[nEntry-1] = 0;
  I[N] = nEntry+1;
  
  //ConjugateGradientCuda * linSolver = new ConjugateGradientCuda();
  //linSolver->initCuda(N,nz, &(I[0]), &(J[0]));
  //linSolver->solve(&(val[0]), &(x[0]), &(rhs[0]));
  //linSolver->clearCuda();

  for(unsigned int ii = 0;ii<rhs.size();ii++){
  //  std::cout<<x[ii]<<"\n";
  }
}

void stiffnessTest()
{
  float h = 0.001f;
  ElementRegGrid * grid = new ElementRegGrid(1, 2, 1);

  //StrainEneNeo* ene=new StrainEneNeo();
  StrainEne* ene = new StrainLin();
  ene->param[1]= 1000;
  ene->param[0] = 10000;
  MaterialQuad * material = new MaterialQuad(ene);
  grid->m.push_back(material);
    grid->x[1][2] +=0.1;
    grid->x[3][1] +=0.2;
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
  //StrainEneNeo ene;
  StrainLin ene;
  ene.param[0] = 1000;
  ene.param[1] = 10000;
  MaterialQuad material(&ene);
  em->m.push_back(&material);
  
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