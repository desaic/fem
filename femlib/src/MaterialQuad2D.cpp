#include "MaterialQuad2D.h"
#include "Element2D.h"
#include "ElementMesh2D.h"
#include "MatrixX.hpp"
#include "Quadrature2D.h"
#include "StrainEne2D.h"
#include <assert.h>

#include <Eigen/Dense>

MaterialQuad2D::MaterialQuad2D(StrainEne2D * ene, Quadrature2D * _q ):q(_q)
{
  if(q==0){
    q = &(Quadrature2D::Gauss2);
  }
  e.resize(q->x.size(), ene);
}

MaterialQuad2D::MaterialQuad2D(const std::vector<StrainEne2D *> & ene, Quadrature2D * _q )
{
  if(q==0){
    q = &(Quadrature2D::Gauss2);
  }
  e=ene;
}

void MaterialQuad2D::init(ElementMesh2D * m)
{
  Element2D * e = m->e[0];
  gradN.resize(e->nV());
  for (int ii = 0; ii < e->nV(); ii++){
    gradN[ii].resize(q->x.size());
    for (int qq = 0; qq < q->x.size(); qq++){
      gradN[ii][qq] = e->shapeFunGrad(ii, q->x[qq], m->X);
    }
  }
}

cfgScalar MaterialQuad2D::getEnergy(Element2D* ele, ElementMesh2D * mesh)
{
  cfgScalar energy = 0;
  for(int ii = 0; ii<q->x.size();ii++){
    Matrix2S F = ele->defGrad(q->x[ii], mesh->X, mesh->x);
    energy += q->w[ii] * e[ii]->getEnergy(F);
  }
  return ele->getVol(mesh->X) * energy;
}

std::vector<MatrixXS> MaterialQuad2D::getElasticityTensors()
{
  std::vector<MatrixXS> elasticityTensors;
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    MatrixXS elasticityTensor = e[ii]->getElasticityTensor();
    elasticityTensors.push_back(elasticityTensor);
  }
  return elasticityTensors;
}

std::vector<Matrix2S> MaterialQuad2D::getStressTensors(Element2D* ele, ElementMesh2D * mesh)
{
  std::vector<Matrix2S> stressTensors;
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    Matrix2S F = ele->defGrad(q->x[ii],mesh->X, mesh->x);
    Matrix2S stressTensor = e[ii]->getStressTensor(F);
    stressTensors.push_back(stressTensor);
  }
  return stressTensors;
}

std::vector<Matrix2S> MaterialQuad2D::getStrainTensors(Element2D* ele, ElementMesh2D * mesh)
{
  return getStrainTensors(ele, mesh, mesh->x);
}

std::vector<Matrix2S> MaterialQuad2D::getStrainTensors(Element2D* ele, ElementMesh2D * mesh, const std::vector<Vector2S> &ix)
{
  assert(ix.size()==mesh->x.size());

  std::vector<Matrix2S> strainTensors;
  for(unsigned int ii = 0; ii<q->x.size(); ii++)
  {
    Matrix2S F = ele->defGrad(q->x[ii],mesh->X, ix);
    Matrix2S strainTensor = e[ii]->getStrainTensor(F);
    strainTensors.push_back(strainTensor);
  }
  return strainTensors;
}

std::vector<Vector2S> MaterialQuad2D::getForce(Element2D* ele, ElementMesh2D * mesh)
{
  std::vector<Vector2S> f(ele->nV(), Vector2S::Zero());
  std::vector<Matrix2S> P(q->w.size());
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    Matrix2S F = ele->defGrad(q->x[ii],mesh->X, mesh->x);
    P[ii] = e[ii]->getPK1(F);
  }
  
  cfgScalar vol = ele->getVol(mesh->X);
  for(unsigned int jj = 0; jj<q->x.size(); jj++){
    for(int ii = 0; ii<ele->nV(); ii++){
      f[ii] -= vol * q->w[jj] * (P[jj]*gradN[ii][jj]);
    }
  }
  return f;
}

MatrixXS MaterialQuad2D::getStiffness(Element2D* ele, ElementMesh2D * mesh)
{
  int ndof = 2* ele->nV();
  MatrixXS K = MatrixXS::Zero(ndof, ndof);

  for(unsigned int ii = 0; ii<q->x.size();ii++)
  {
    K += q->w[ii] * stiffness(ii, this, ele, mesh);
  }
  cfgScalar vol = ele->getVol(mesh->X);
  K *= vol;
  MatrixXS Kret(ndof, ndof);
  copyMat(K, Kret, ndof, ndof);
  return Kret;
}

MatrixXS MaterialQuad2D::stiffness(int qi, const MaterialQuad2D * mat, Element2D* ele, ElementMesh2D * mesh)
{
  int nquad = (int)mat->q->x.size();
  int ndof = 2*ele->nV();
  Vector2S p = mat->q->x[qi];

  MatrixXS K = MatrixXS::Zero(ndof, ndof);
  Matrix2S F = ele->defGrad(p,mesh->X,mesh->x);
  
  for(int ii = 0;ii<4;ii++){
    for(int jj = 0;jj<2;jj++){
      Matrix2S dF = Matrix2S::Zero();
      dF.row(jj) = mat->gradN[ii][qi];
      //dF.setRow(jj, mat->gradN[ii][qi]);
      Matrix2S dP = mat->e[ii]->getdPdx(F,dF);
      for(int vv = 0;vv<4;vv++){
        Vector2S dfdxi = dP*mat->gradN[vv][qi];
        int col = 2*ii+jj;
        for(int kk = 0;kk<2;kk++){
          K(2*vv+kk, col) = dfdxi[kk];
        }
      }
    }
  }
  return K;
}
