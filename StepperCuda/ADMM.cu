#include "ADMM.h"
#include "HexEle.h"

#include "ElementMesh.hpp"
#include "Element.hpp"
#include "vecmath.h"
#include "Util.hpp"
#include <time.h>  
#include <fstream>
ADMMStepper::ADMMStepper():NSteps(100),nThread(256),ro(100),tol(0.0001f),
  maxDist(0.05f),roMult(1.5f),e(0), outname("out.txt")
{}

ADMMStepper::~ADMMStepper()
{
  ///@TODO delete arrays and cuda free
}

float3 vec2float(const Vector3f & v)
{
  return make_float3(v[0],v[1],v[2]);
}

Vector3f float2vec(const float3 & v)
{
  return Vector3f(v.x, v.y, v.z);
}

void ADMMStepper::setZdev(const std::vector<Vector3f> & vec)
{
  int nEle = (int)e->e.size();
  for(int ii = 0;ii<nEle;ii++){
    ADMMInfo * admm = &(hostadmm[ii]);
    for(int jj = 0;jj<NVERT;jj++){
      int vIdx = e->e[ii]->at(jj);
      admm->Z[jj] = vec2float(vec[vIdx]);
    }
  }
  cudaMemcpy(devadmm, hostadmm, nEle*sizeof(ADMMInfo), cudaMemcpyHostToDevice);
}

float ADMMStepper::getEnergy()
{
  int nEle = (int)e->e.size();
  int nBlock = nEle/nThread + ((nEle % nThread)!=0);
  GetEnergy<<<nBlock, nThread>>>(devXX, devadmm, Edev);
  cudaMemcpy(Ehost, Edev, nEle*sizeof(float), cudaMemcpyDeviceToHost);
  float ene = 0;
  for(int ii= 0;ii<nEle;ii++){
    ene += Ehost[ii];
  }
  for(unsigned int ii = 0;ii<e->fe.size();ii++){
    ene -= Vector3f::dot(e->fe[ii], Z[ii]);
  }
  return ene;
}

//get forces for the mesh excluding ADMM terms
void ADMMStepper::getForce(std::vector<Vector3f> & ff)
{
  int nEle = (int)e->e.size();
  int nBlock = nEle/nThread + ((nEle % nThread)!=0);
  
  GetIntForce<<<nBlock, nThread>>>(devXX, devadmm, fdev);
  cudaMemcpy(fhost, fdev, nEle*NVERT*sizeof(float3), cudaMemcpyDeviceToHost);
  //cumulate internal forces
  ff.assign(ff.size(), Vector3f::ZERO);
  for(int ii= 0;ii<nEle;ii++){
    for(int jj = 0;jj<NVERT;jj++){
      int vidx = e->e[ii]->at(jj);
      float3 ftmp = fhost[NVERT*ii + jj];
      Vector3f fint(ftmp.x,ftmp.y,ftmp.z);
      ff[vidx] += fint;
    }
  }
  for(int ii = 0;ii<e->fe.size();ii++){
    ff[ii] += e->fe[ii];
  }
  for(unsigned int ii = 0;ii<e->fixed.size();ii++){
    if (e->fixed[ii]){
      ff[ii] = Vector3f::ZERO;
    }
  }
}

void ADMMStepper::stepGrad()
{
  int nEle = (int)e->e.size();
  
  //change in Z
  std::vector<Vector3f> dZ(Z.size());
  //internal + external forces
  std::vector<Vector3f> ff(Z.size());
  float hh = 0.001f;
  for(int iter = 0;iter<NSteps;iter++){
    std::cout<<iter<<"\n";
    setZdev(Z);
    float ene = getEnergy();
    getForce(ff);
    float mag = 0;
    for(unsigned int ii = 0;ii<ff.size(); ii++){
      for(int jj = 0;jj<3;jj++){
        mag += std::abs(ff[ii][jj]);
      }
    }
    mag = mag/ff.size();
    if(mag<tol){
      break;
    }
    float ene1=ene;
    while(1){
      dZ = hh * ff;
      Z = e->x + dZ;
      setZdev(Z);
      ene1 = getEnergy();
      if(ene1<ene){
        hh=hh*1.5f;
        break;
      }else{
        hh=hh/2;
        if(hh<1e-15f){
          break;
        }
      }
    }
    std::cout<<hh<<" "<<ene1<<"\n";
    
    e->x = Z;
  }
}

void
ADMMStepper::step()
{
  int nEle = (int)e->e.size();
  int nBlock = nEle/nThread + ((nEle % nThread)!=0);
  
  clock_t tt0, tt;
  tt0 = clock();
  std::ofstream out(outname);

  //change in Z
  std::vector<Vector3f> dZ(Z.size());
  //internal + external forces
  std::vector<Vector3f> ff(Z.size());

  for(int iter = 0;iter<NSteps;iter++){
    std::cout<<iter<<"\n";
    tt = clock();
    out<<(tt-tt0)/(CLOCKS_PER_SEC/(float)1000)<<" ";
    setZdev(Z);
    float ene = getEnergy();
    out<<ene<<"\n";
    getForce(ff);
    float eleSize = e->eleSize();
    float maxRo = 0;
    float maxDiff=0;
    //adjust ro
    for(int ii = 0;ii<nEle;ii++){
      ADMMInfo * admm = &(hostadmm[ii]);
      for(int jj = 0;jj<NVERT;jj++){
        int vidx = e->e[ii]->at(jj);
        Vector3f zz = Z[vidx];
        Vector3f xx = float2vec(hostxx[ii*NVERT+jj]);
        float diff = (xx-zz).abs();
        if(diff>maxDiff){
          maxDiff = diff;
        }
        if(admm->ro > maxRo){
          maxRo = admm->ro;
        }
        if( diff > maxDist*eleSize){
          admm->ro *= roMult;
          break;
        }
      }
    }
    std::cout<<"maxro "<<maxRo<<"\n";
    std::cout<<"maxdiff "<<maxDiff<<"\n";
    for(int ii = 0;ii<nEle;ii++){
      //copy Z and y to device
      ADMMInfo * admm = &(hostadmm[ii]);
      for(int jj = 0;jj<NVERT;jj++){
        int vIdx = e->e[ii]->at(jj);
        admm->Z[jj] = vec2float(Z[vIdx]);
        admm->y[jj] = vec2float(l[ii][jj]);
      }
    }
    cudaMemcpy(devadmm, hostadmm, nEle*sizeof(ADMMInfo), cudaMemcpyHostToDevice);
    admmMinEleDup<<<nBlock, nThread>>>(devXX, devxx, devadmm);
    cudaMemcpy(hostxx, devxx, nEle*NVERT*sizeof(float3), cudaMemcpyDeviceToHost);
       
    //update z closed form
    for(size_t ii = 0;ii<e->X.size();ii++){
		  Z[ii] = Vector3f::ZERO;
    }
    //add per element variables and multipliers
    for(int ii = 0;ii<nEle;ii++){
      Element * ele = e->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        Vector3f xx = float2vec(hostxx[ii*NVERT+jj]);
        Z[ele->at(jj)] += xx + l[ii][jj]/ro;
      }
    }

    //add force variables
    for(size_t ii = 0;ii<e->fe.size();ii++){
      Z[ii] += (1.0f/ro)*e->fe[ii];
    }

    //divide
    for(size_t ii = 0;ii<e->x.size();ii++){
      Z[ii] /= (float)N[ii];
    }
    
    //fix constrained vertices
    for(auto ii = 0;ii<e->fixed.size();ii++){
      if (e->fixed[ii]){
        Z[ii] = e->x[ii];
      }
    }

    //update multiplier for elements
    for(size_t ii = 0;ii<nEle;ii++){
      Element * ele = e->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
         Vector3f xx = float2vec(hostxx[ii*NVERT+jj]);
        l[ii][jj] += ro*(xx - Z[ele->at(jj)]);
      }
    }

    float hh = 1;
    float mag = 0;
    for(unsigned int ii = 0;ii<ff.size(); ii++){
      for(int jj = 0;jj<3;jj++){
        mag += std::abs(ff[ii][jj]);
      }
    }
    mag = mag/ff.size();
    if(mag<tol){
      break;
    }
    
    float ene1=ene;
    for(unsigned int ii =0;ii<Z.size();ii++){
      dZ[ii] = Z[ii] - e->x[ii];
    }
    while(1){
      Z = e->x + hh*dZ;
      setZdev(Z);
      ene1 = getEnergy();
      if(ene1<ene){
        break;
      }else{
        hh=hh/2;
        if(hh<1e-15){
          break;
        }
      }
    }
    std::cout<<hh<<" "<<ene1<<"\n";
    e->x = Z;
  }
}

void
ADMMStepper::initVar(ElementMesh * _e)
{
  e=_e;
  int nEle = (int)e->e.size();
  l.resize(e->e.size());
  N.assign(e->X.size(),0);
  for(size_t ii= 0; ii<e->e.size();ii++){
    Element * ele = e->e[ii];
    l[ii].resize(ele->nV(), Vector3f::ZERO);
    for(int jj = 0;jj<ele->nV();jj++){
      N[ele->at(jj)] ++;
    }
  }

  Z = e->x;
  int xsize = NVERT * nEle;
  hostXX = new float3[NVERT];
  hostxx = new float3 [xsize];
  Element * ele = e->e[0];
  
  for(int ii = 0;ii<NVERT;ii++){
    hostXX[ii].x= e->X[ele->at(ii)][0];
	  hostXX[ii].y= e->X[ele->at(ii)][1];
	  hostXX[ii].z= e->X[ele->at(ii)][2];
  }
  hostadmm = new ADMMInfo[nEle];
  cudaMalloc((void**)&devXX,        NVERT*sizeof(float3));
  cudaMalloc((void**)&devxx,   nEle*NVERT*sizeof(float3));
  cudaMalloc((void**)&devadmm, nEle*sizeof(ADMMInfo));
  cudaMemcpy(devXX,hostXX, NVERT*sizeof(float3),cudaMemcpyHostToDevice);

  int xIdx = 0;
  for(size_t ii = 0;ii<e->e.size();ii++){
    Element * ele = e->e[ii];
    hostadmm[ii].ro = ro;
    for(int jj = 0;jj<ele->nV();jj++){
      hostxx[xIdx] = vec2float(e->x[ele->at(jj)]);
      xIdx++;
    }
  }
  cudaMemcpy(devxx, hostxx, xsize * sizeof(float3), cudaMemcpyHostToDevice);
  initHexEle();

  Ehost = new float[nEle];
  cudaMalloc(&Edev,nEle * sizeof(float));
  fhost = new float3[nEle*NVERT];
  cudaMalloc(&fdev, nEle*NVERT * sizeof(float3));

}
