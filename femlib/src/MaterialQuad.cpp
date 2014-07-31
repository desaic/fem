#include "MaterialQuad.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "MatrixXd.hpp"
#include "Quadrature.hpp"
#include "StrainEne.hpp"

#include <Eigen/Dense>

///@brief helper for computing stiffness contribution of one quadrature point
Eigen::MatrixXf stiffness(int qi, const MaterialQuad* mat, Element* ele, ElementMesh * mesh);

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
  int ndof = 3* ele->nV();
  Eigen::MatrixXf K = Eigen::MatrixXf::Zero(ndof, ndof);

  for(unsigned int ii = 0; ii<q->x.size();ii++){
    K += q->w[ii] * stiffness(ii, this, ele, mesh);
  }
  float vol = ele->getVol(mesh->X);
  K *= vol;
  MatrixXd Kret(ndof, ndof);
  copyMat(K, Kret, ndof, ndof);
  return Kret;
}

Eigen::MatrixXf getStiffness(int qi, const MaterialQuad * mat, Element* ele, ElementMesh * mesh)
{
  int nquad = (int)mat->q->x.size();
  int ndof = 3*ele->nV();
  Vector3f p = mat->q->x[qi];

  std::vector<Vector3f> dN(nquad);
  Eigen::MatrixXf K = Eigen::MatrixXf::Zero(ndof, ndof);
  Matrix3f F = ele->defGrad(p,mesh->X,mesh->x);
  for(int ii = 0;ii<8;ii++){
    dN[ii] = ele->shapeFunGrad(ii,p,mesh->X);
  }
  for(int ii = 0;ii<8;ii++){
    for(int jj = 0;jj<3;jj++){
      Matrix3f dF;
      dF.setRow(jj, dN[ii]);
      Matrix3f dP = mat->e[ii]->getdPdx(F,dF);
      for(int vv = 0;vv<8;vv++){
        Vector3f dfdxi = dP*dN[vv];
        int col = 3*ii+jj;
        for(int kk = 0;kk<3;kk++){
          K(3*vv+kk, col) = dfdxi[kk];
        }
      }
    }
  }
  return K;
}
