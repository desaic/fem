#include "MaterialQuad.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "Quadrature.hpp"
#include "StrainEne.hpp"

#include <Eigen/Dense>
using namespace Eigen;
///@brief helper for computing stiffness contribution of one quadrature point
MatrixXS stiffness(int qi, const MaterialQuad * mat, Element* ele, ElementMesh * mesh);

MaterialQuad::MaterialQuad(StrainEne * ene, Quadrature * _q ):q(_q)
{
  if(q==0){
    q = &(Quadrature::Gauss2);
  }
  e.resize(q->x.size(), ene);
}

MaterialQuad::MaterialQuad(const std::vector<StrainEne *> & ene, Quadrature * _q )
{
  if(q==0){
    q = &(Quadrature::Gauss2);
  }
  e=ene;
}

void MaterialQuad::init(ElementMesh * m)
{
  Element * e = m->e[0];
  gradN.resize(e->nV());
  for (int ii = 0; ii < e->nV(); ii++){
    gradN[ii].resize(q->x.size());
    for (int qq = 0; qq < q->x.size(); qq++){
      gradN[ii][qq] = e->shapeFunGrad(ii, q->x[qq], m->X);
    }
  }
}

cfgScalar MaterialQuad::getEnergy(Element* ele, ElementMesh * mesh)
{
  cfgScalar energy = 0;
  for(int ii = 0; ii<q->x.size();ii++){
    Matrix3S F = ele->defGrad(q->x[ii], mesh->X, mesh->x);
    energy += q->w[ii] * e[ii]->getEnergy(F);
  }
  return ele->getVol(mesh->X) * energy;
}

std::vector<Vector3S> MaterialQuad::getForce(Element* ele, ElementMesh * mesh)
{
  std::vector<Vector3S> f(ele->nV(), Vector3S::Zero());
  std::vector<Matrix3S> P(q->w.size());
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    Matrix3S F = ele->defGrad(q->x[ii],mesh->X, mesh->x);
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

MatrixXS MaterialQuad::getStiffness(Element* ele, ElementMesh * mesh)
{
  int ndof = 3* ele->nV();
  MatrixXS K = MatrixXS::Zero(ndof, ndof);

  for(unsigned int ii = 0; ii<q->x.size();ii++){
    K += q->w[ii] * stiffness(ii, this, ele, mesh);
  }
  cfgScalar vol = ele->getVol(mesh->X);
  K *= vol;
  return K;
}

MatrixXS stiffness(int qi, const MaterialQuad * mat, Element* ele, ElementMesh * mesh)
{
  int nquad = (int)mat->q->x.size();
  int ndof = 3*ele->nV();
  Vector3S p = mat->q->x[qi];

  MatrixXS K = MatrixXS::Zero(ndof, ndof);
  Matrix3S F = ele->defGrad(p,mesh->X,mesh->x);
  
  for(int ii = 0;ii<8;ii++){
    for(int jj = 0;jj<3;jj++){
      Matrix3S dF=Matrix3S::Zero();
      dF.row(jj)=mat->gradN[ii][qi];
      Matrix3S dP = mat->e[ii]->getdPdx(F,dF);
      for(int vv = 0;vv<8;vv++){
        Vector3S dfdxi = dP*mat->gradN[vv][qi];
        int col = 3*ii+jj;
        for(int kk = 0;kk<3;kk++){
          K(3*vv+kk, col) = dfdxi[kk];
        }
      }
    }
  }
  return K;
}

std::vector<MatrixXS> MaterialQuad::getElasticityTensors()
{
  std::vector<MatrixXS> elasticityTensors;
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    MatrixXS elasticityTensor = e[ii]->EMatrix();
    elasticityTensors.push_back(elasticityTensor);
  }
  return elasticityTensors;
}

std::vector<Matrix3S> MaterialQuad::getStrainTensors(Element* ele, ElementMesh * mesh, const std::vector<Vector3S> &ix)
{
  assert(ix.size()==mesh->x.size());

  std::vector<Matrix3S> strainTensors;
  for(unsigned int ii = 0; ii<q->x.size(); ii++)
  {
    Matrix3S F = ele->defGrad(q->x[ii],mesh->X, ix);
    Matrix3S strainTensor = e[ii]->getStrainTensor(F);
    strainTensors.push_back(strainTensor);
  }
  return strainTensors;
}


