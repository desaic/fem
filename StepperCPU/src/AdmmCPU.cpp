#include "AdmmCPU.hpp"
#include "Element.hpp"
#include "ElementHex.hpp"
#include "ElementMesh.hpp"
#include "MatrixXd.hpp"
#include "LinSolve.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"

#include <time.h>
#include <fstream>

float AdmmCPU::getEnergy(ElementMesh * eMesh, int eIdx)
{
  Element * e = eMesh->e[eIdx];
  float E = eMesh->getEnergy(eIdx);
  for(int ii = 0; ii<e->nV(); ii++){
    Vector3f diff = eMesh->x[e->at(ii)]-Z[e->at(ii)];
    E += 0.5f * ro[eIdx] * Vector3f::dot(diff,diff);
    E += Vector3f::dot(y[eIdx][ii], diff);
  }
  return E;
}

std::vector<Vector3f>
  AdmmCPU::getForces(ElementMesh * eMesh, int eIdx)
{
  Element * ele = eMesh->e[eIdx];
  std::vector<Vector3f> ff = eMesh->getForce(eIdx);
  for(unsigned int ii = 0;ii<ff.size();ii++){
    Vector3f diff = eMesh->x[ele->at(ii)] - Z[ele->at(ii)];
    ff[ii] -= ro[eIdx] * diff;
    ff[ii] -= y[eIdx][ii];
  }
  return ff;
}

MatrixXd
  AdmmCPU::stiffness(ElementMesh *mesh, int eIdx)
{
  MatrixXd K = mesh->getStiffness(eIdx);
  for(int ii = 0;ii<K.mm;ii++){
    K(ii,ii) += ro[eIdx];
  }
  return K;
}

void AdmmCPU::minimizeElement(ElementMesh * m, Element * ele,
                              int eIdx)
{
  float E = 0;
  float h = 1;
  int NSteps = 100;
  int ndof = 3*ele->nV();
  
  for(int iter = 0; iter<NSteps; iter++){
    std::vector<Vector3f> force = getForces(m,eIdx);
    MatrixXd K = stiffness(m,eIdx);
    K.print(std::cout);
    for(int ii = 0; ii<ele->nV(); ii++){
      int vidx = ele->at(ii);
      if(m->fixed[vidx]){
        force[ii] = Vector3f(0,0,0);
      }
      for(int jj = 0;jj<3;jj++){
        int row = 3*ii + jj;
        K(row,row) += 100;
        bb[ row ] = force[ii][jj];
        if(m->fixed[ii]){
          for(int kk = 0;kk<ndof;kk++){
            K(row,kk) = 0;
            K(row,row) = 1;
          }
        }
      }
    }
    linSolve(K,bb);

    for(int ii = 0; ii<ele->nV(); ii++){
      for(int jj = 0;jj<3;jj++){
        force[ii][jj] = (float)bb[3*ii+jj];
      }
    }

    float totalMag = 0;
    for(unsigned int ii = 0;ii<force.size();ii++){
      int vidx = ele->at(ii);
      for(int jj =0 ; jj<3; jj++){
        totalMag += std::abs(force[ii][jj]);
      }
    }
    if(totalMag<xtol){
      return ;
    }

    float E = getEnergy(m,eIdx);
    
    //line search
    std::vector<Vector3f> x0 = m->x;
    float E1;
    while(1){
      m->x=x0;
      addmul(m->x, h, force);
      E1 = getEnergy(m,eIdx);
      std::cout<<h<<" "<<E1<<"\n";
      if(E1>E || fem_error){
        fem_error = 0;
        h = 0.5f* h;
      }else{
        h=1.1f*h;
        break;
      }
    }
  }
}

void AdmmCPU::initVar(ElementMesh *e)
{
  if(bb!=0){
    delete[] bb;
  }
  bb = new double[3*e->x.size()];
  u.resize(e->e.size());
  y.resize(u.size());
  ro.resize(e->e.size(),ro0);
  for(size_t ii= 0; ii<e->e.size();ii++){
    Element * ele = e->e[ii];
    u[ii].resize(ele->nV());
    y[ii].resize(u[ii].size());
    for(int jj = 0;jj<ele->nV();jj++){
      u[ii][jj] = e->X[ele->at(jj)];
    }
  }
  Z = e->x;
}

void AdmmCPU::step(ElementMesh * m)
{
  initVar(m);
  
  std::ofstream out("converge.txt");
  clock_t tt,tt0 ;
  tt0 = clock();

  std::vector<Vector3f> x0 = m->x;

  bool pause = true;
  m->u = &u;
  float prevE = m->getEnergy();
  
  for(int iter = 0;iter<nSteps;iter++){
    float maxRo = 0;
    float maxDiff=0;
    //adjust ro
    for(unsigned int ii = 0;ii<m->e.size();ii++){
      Element * ele = m->e[ii];
      float eleSize = m->X[ele->at(1)][2] - m->X[ele->at(0)][2];
      for(int jj = 0; jj<ele->nV(); jj++){
        int vidx = ele->at(jj);
        Vector3f zz = Z[vidx];
        Vector3f xx = u[ii][jj];
        float diff = (xx-zz).abs();
        if(diff>maxDiff){
          maxDiff = diff;
        }
        if(ro[ii] > maxRo){
          maxRo = ro[ii];
        }
        if( diff > maxDist*eleSize){
          ro[ii] *= roMult;
          break;
        }
      }
    }
    std::cout<<"maxro "<<maxRo<<"\n";
    std::cout<<"maxdiff "<<maxDiff<<"\n";

    //Z in the previous iteraiton.
    std::vector<Vector3f> Zk_1=Z;

    //update u locally
    for(unsigned int ii = 0;ii<m->e.size();ii++){
      Element * ele = m->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        m->x[ele->at(jj)] = u[ii][jj];
      }
      minimizeElement(m,ele, ii);

      for(int jj = 0;jj<ele->nV();jj++){
        u[ii][jj] = m->x[ele->at(jj)];
      }
    }

    //update z closed form
    for(unsigned int ii = 0;ii<m->X.size();ii++){
      Z[ii] = Vector3f(0,0,0);
    }

    std::vector<float> roSum(Z.size(),0.0f);

    //add per element variables and multipliers
    for(unsigned int ii = 0;ii<u.size();ii++){
      Element * ele = m->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        int vIdx = ele->at(jj);
        Z[vIdx] += ro[ii] * u[ii][jj]+y[ii][jj];
        roSum[vIdx] += ro[ii];
      }
    }

    for(unsigned int ii = 0; ii< m->fe.size();ii++){
      Vector3f force = m->fe[ii];
      Z[ii] += force;
    }
    
    for(size_t ii = 0;ii<m->x.size();ii++){
      Z[ii] /= roSum[ii];
    }

    //fix constraints
    for(unsigned int ii = 0;ii<m->fixed.size();ii++){
      if(m->fixed[ii]){
        Z[ii] = x0[ii];
      }
    }

    //update multiplier for elements
    for(unsigned int ii = 0;ii<u.size();ii++){
      Element * ele = m->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        y[ii][jj] +=  ro[ii] * (u[ii][jj] - Z[ele->at(jj)]);
      }
    }
    float E = m->getEnergy();
    std::cout<<"Energy in iteration "<<iter<<": "<<E<<"\n";
    
    float ene1=E;
    for(unsigned int ii =0;ii<Z.size();ii++){
      Z[ii] = Z[ii] - Zk_1[ii];
    }
    float hh = 1.0f;
    while(1){
      m->x=Zk_1;
      addmul(m->x,hh,Z);
      ene1 = m->getEnergy();
      if(ene1<E && fem_error==0){
        Z=m->x;
        break;
      }else{
        hh=hh/2;
        if(hh<1e-15){
          m->x = Zk_1;
          Z=Zk_1;
          break;
        }
      }
    }
    std::cout<<"hh "<<hh<<"\n";    
    if(prevE-ene1 < tol){
      break;
    }
    prevE = E;

    tt = clock();
    out<<(tt-tt0)/(CLOCKS_PER_SEC/1000.0);

  }
}

AdmmCPU::~AdmmCPU()
{
  if(bb!=0){
    delete bb;
  }
}

AdmmCPU::AdmmCPU():bb(0),maxDist(0.05f),ro0(1000.0f),
  roMult(1.5f),tol(0.001f),xtol(0.001f)
{
}
