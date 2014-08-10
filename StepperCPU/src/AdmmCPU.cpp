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
    Vector3f diff = eMesh->x[e->at(ii)]-Z[e->at(ii)] + y[eIdx][ii];
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
    Vector3f diff = eMesh->x[ele->at(ii)] -Z[ele->at(ii)]+y[eIdx][ii];
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

std::vector<Vector3f>
  getEleX(ElementMesh * mesh, Element * ele)
{
  std::vector<Vector3f> xx(ele->nV());
  for(int ii = 0; ii<ele->nV(); ii++){
    xx[ii] = mesh->x[ele->at(ii)];
  }
  return xx;
}

std::vector<Vector3f>
  setEleX(ElementMesh * mesh, Element * ele, const std::vector<Vector3f> & xx)
{
  for(int ii = 0;ii<ele->nV();ii++){
    mesh->x[ele->at(ii)] = xx[ii];
  }
}

void AdmmCPU::minimizeElement(ElementMesh * m, Element * ele,
                              int eIdx)
{
  float E = 0;
  int ii = 0;
  float h = 1;
  int NSteps = 100;
  int ndof = 3*ele->nV();
  double * bb = new double[ndof];

  for(ii = 0;ii<NSteps;ii++){
    std::vector<Vector3f> force = m->getForce();
    float E = m->getEnergy();
    float totalMag = 0;
    for(unsigned int ii = 0;ii<force.size();ii++){
      totalMag += force[ii].absSquared();  
    }
    if(totalMag<xtol){
      return ;
    }

    MatrixXd K = m->getStiffness();

    int ndof = 3*(int)m->x.size();
    for(unsigned int ii = 0;ii<m->x.size(); ii++){
      for(int jj = 0;jj<3;jj++){
        int row = 3*ii + jj;
        K(row,row) += 100;
        if(m->fixed[ii]){
          bb[ row ] = 0;
          for(int kk = 0;kk<ndof;kk++){
            K(row,kk) = 0;
            K(row,row) = 1;
          }
        }else{
          bb[ row ] = force[ii][jj];
        }
      }
    }
    linSolve(K,bb);

    for(unsigned int ii = 0;ii<m->x.size(); ii++){
      for(int jj = 0;jj<3;jj++){
        force[ii][jj] = (float)bb[3*ii+jj];
      }
    }

    //line search
    std::vector<Vector3f> x0 = m->x;
    float E1;
    while(1){
      m->x=x0;
      addmul(m->x, h, force);
      E1 = m->getEnergy();

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
  u.resize(e->e.size());
  y.resize(u.size());
  N.assign(e->X.size(),0);
  for(size_t ii= 0; ii<e->e.size();ii++){
    Element * ele = e->e[ii];
    u[ii].resize(ele->nV());
    y[ii].resize(u[ii].size());
    for(int jj = 0;jj<ele->nV();jj++){
      u[ii][jj] = e->X[ele->at(jj)];
      N[ele->at(jj)] ++;
    }
  }
  e->x = e->X;
  Z = e->x;
}

void AdmmCPU::step(ElementMesh * e)
{
  initVar(e);
  float E = e->getEnergy();
  float prevE = E;

  std::ofstream out("converge.txt");
  clock_t tt,tt0 ;
  tt0 = clock();

  std::vector<Vector3f> Z_k = e->x;

  bool pause = true;
  e->u = &u;
  for(int iter = 0;iter<nSteps;iter++){
    //objective. Augmentedd lagrangian
    float obj = 0;
    //update u locally
    for(unsigned int ii = 0;ii<e->e.size();ii++){
      Element * ele = e->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        e->x[ele->at(jj)] = u[ii][jj]; //Z[ele->At(jj)]
      }
      minimizeElement(e,ele, ii);

      for(int jj = 0;jj<ele->nV();jj++){
        u[ii][jj] = e->x[ele->at(jj)];
      }
    }

    //update z closed form
    for(unsigned int ii = 0;ii<e->X.size();ii++){
      Z[ii] = Vector3f(0,0,0);
    }

    //add per element variables and multipliers
    for(unsigned int ii = 0;ii<u.size();ii++){
      Element * ele = e->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        Z[ele->at(jj)] += u[ii][jj]+y[ii][jj];
      }
    }

    //fix constraints
    for(unsigned int ii = 0;ii<e->fixed.size();ii++){
    }

    //divide
    for(size_t ii = 0;ii<e->x.size();ii++){
      Z[ii] /= (float)N[ii];
    }

    //update multiplier for elements
    for(unsigned int ii = 0;ii<u.size();ii++){
      Element * ele = e->e[ii];
      for(int jj = 0;jj<ele->nV();jj++){
        y[ii][jj] +=  u[ii][jj] - Z[ele->at(jj)];
      }
    }
    Z_k = Z;
    e->x = Z;
    E = e->getEnergy();
    std::cout<<"Energy in iteration "<<iter<<": "<<E<<"\n";
    
    if(std::abs(prevE-E) < tol){
      break;
    }
    prevE = E;

    tt = clock();
    out<<(tt-tt0)/(CLOCKS_PER_SEC/1000.0);
    float ene = e->getEnergy();
    out<<" "<<ene<<"\n";
  }
}

AdmmCPU::AdmmCPU():tol(0.001f),xtol(0.001f)
{
}
