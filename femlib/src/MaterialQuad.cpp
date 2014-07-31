#include "MaterialQuad.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "MatrixXd.hpp"
#include "Quadrature.hpp"
#include "StrainEne.hpp"

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

float MaterialQuad::getEnergy(Element* ele, ElementMesh * mesh)
{
  float energy = 0;
  for(int ii = 0; ii<q->x.size();ii++){
    Matrix3f F = ele->defGrad(q->x[ii], mesh->X, mesh->x);
    energy += q->w[ii] * e[ii]->getEnergy(F);
  }
  return ele->getVol(mesh->X) * energy;
}

std::vector<Vector3f> MaterialQuad::getForce(Element* ele, ElementMesh * mesh)
{
  std::vector<Vector3f> f(ele->nV(), Vector3f::ZERO);
  std::vector<Matrix3f> P(q->w.size());
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    Matrix3f F = ele->defGrad(q->x[ii],mesh->X, mesh->x);
    P[ii] = e[ii]->getPK1(F);
  }
  
  float vol = ele->getVol(mesh->X);
  for(unsigned int jj = 0; jj<q->x.size(); jj++){
    for(int ii = 0; ii<ele->nV(); ii++){
      Vector3f gradN = ele->shapeFunGrad(ii, q->x[jj], mesh->X);
      f[ii] -= vol * q->w[jj] * (P[jj]*gradN);
    }
  }
  return f;
}

MatrixXd MaterialQuad::getStiffness(Element* ele, ElementMesh * mesh)
{
 
}