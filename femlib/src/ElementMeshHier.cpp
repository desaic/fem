#include "ElementHier.hpp"
#include "ElementMeshHier.hpp"
#include "ElementMesh.hpp"
#include "ElementRegGrid.hpp"
#include "Element.hpp"
#include "MaterialQuad.hpp"
#include "Quadrature.hpp"
#include "StrainEne.hpp"
#include <iostream>
const float eps = 1e-6f;
///@brief construct nLevel additional levels given m[0].
void ElementMeshHier::buildHier(int nLevel)
{
  if (m.size() == 0){
    std::cout << "Error ElementMeshHier coarsen. Need fine mesh.\n";
    return;
  }
  if (nLevel < 1){
    return;
  }
  replaceElementHex((ElementRegGrid*)m[0]);
  for (int ii = 1; ii <= nLevel; ii++){
    //fine mesh
    ElementRegGrid * fm = (ElementRegGrid*)(m[ii - 1]);
    //coarse mesh
    ElementRegGrid * cm = coarsen(fm);
    m.push_back((ElementMesh*)cm);
    std::vector<bool> visited(fm->x.size(), false);
    for (unsigned int ei = 0; ei < fm->e.size(); ei++){
      //fine element
      ElementHier * ef = (ElementHier*)(fm->e[ei]);
      int parent = ef->parent;
      //coarse element
      ElementHier * ec = (ElementHier*)(cm->e[parent]);
      ec->children.push_back(ei);
      //coarsen forces and constraints
      for (int vi = 0; vi < 8; vi++){
        //fine vertex index
        int fv = ef->at(vi);
        if (visited[fv]){
          continue;
        }
        visited[fv] = true;
        Vector3f nat = ef->Xn[vi];
        std::vector<float> N = ec->shapeFun(nat);
        for (unsigned int jj = 0; jj< N.size(); jj++){
          //coarse vertex
          int cv = ec->at(jj);
          cm->fe[cv] += N[jj] * fm->fe[fv];
          if (fm->fixed[fv] && N[jj] > eps){
            cm->fixed[cv] = true;
          }
        }
      }
    }
  }

  //generate list of quadrature points
  ElementRegGrid * fm = (ElementRegGrid*)(m[0]);
  for (unsigned int ii = 0; ii < fm->e.size(); ii++){
    MaterialQuad * mat = (MaterialQuad * )(fm->m[fm->me[ii]]);
    ElementHier* fe = (ElementHier*)(fm->e[ii]);
    const Quadrature * quadrature = mat->q;
    for (unsigned int qi = 0; qi < quadrature->x.size(); qi++){
      QuadPt q;
      q.F.push_back(Matrix3f::identity());
      q.w.push_back(quadrature->w[qi]);
      q.X.push_back(quadrature->x[qi]);
      q.ei.push_back(ii);
      //coordinate of quadrature point in reference frame
      Vector3f X(0, 0, 0);
      std::vector<float>N = fe->shapeFun(q.X[0]);
      for (unsigned int jj = 0; jj < N.size(); jj++){
        X += N[jj] * fm->X[fe->at(jj)];
      }
      fe->children.push_back((int)quadpt.size());
      int feidx = ii;
      for (unsigned int level = 1; level < m.size(); level++){
        ElementHier* e0 = (ElementHier*)(m[level - 1]->e[feidx]);
        //coarse element
        int cvidx = e0->parent;
        ElementHier * ce = (ElementHier*)(m[level]->e[cvidx]);
        q.F.push_back(Matrix3f::identity());
        //@brief uniformly scale each quadrature point by 1/8;
        q.w.push_back(q.w.back() * 0.125f);
        Vector3f nat = ce->natCoord(X, m[level]->X);
        q.X.push_back(nat);
        q.ei.push_back(cvidx);
        feidx = cvidx;
      }
      
      quadpt.push_back(q);
    }
  }
  //in case input is not at rest pose.
  for (int ii = 0; ii < nLevel; ii++){
    updateDefGrad(ii);
  }
}

void ElementMeshHier::updateDefGrad(int level, int qidx)
{
  int ei = quadpt[qidx].ei[level];
  Element * e = m[level]->e[ei];
  Matrix3f F = e->defGrad(quadpt[qidx].X[level], m[level]->X, m[level]->x);
  quadpt[qidx].F[level] = F;
}

void ElementMeshHier::updateDefGrad(int level)
{
  for (unsigned int ii = 0; ii < quadpt.size(); ii++){
    updateDefGrad(level, ii);
  }
}

ElementMeshHier::ElementMeshHier()
{
}

std::vector<int> ElementMeshHier::getQuadIdx(int level, int eIdx)
{
  //a layer of elements at level k
  std::vector<int> layer;
  layer.push_back(eIdx);
  for (int ii = level - 1; ii >= 0; ii--){
    std::vector<int> layer1;
    for (unsigned int jj = 0; jj < layer.size(); jj++){
      int ei = layer[jj];
      ElementHier * ce = (ElementHier *)(m[ii+1]->e[ei]);
      layer1.insert(layer1.end(), ce->children.begin(), ce->children.end());
    }
    layer = layer1;
  }
  //get list of quadrature points
  std::vector<int> qi;
  for (unsigned int ii = 0; ii < layer.size(); ii++){
    int ei = layer[ii];
    ElementHier * ce = (ElementHier *)(m[0]->e[ei]);
    qi.insert(qi.end(), ce->children.begin(), ce->children.end());
  }
  return qi;
}

float ElementMeshHier::getEnergy()
{
  float E = 0;
  ElementMesh * mesh = m[0];
  for (unsigned int ii = 0; ii < mesh->e.size(); ii++){
    E += getEnergy(0, ii);
  }
  return E;
}

float ElementMeshHier::getEnergy(int level, int eIdx)
{
  float E = 0;
  std::vector<int> quadIdx = getQuadIdx(level, eIdx);
  Element * ele = m[level]->e[eIdx];
  float vol = ele->getVol(m[level]->X);

  for (unsigned int ii = 0; ii < quadIdx.size(); ii++){
    int qi = quadIdx[ii];
    updateDefGrad(level, qi);
    QuadPt * q = &quadpt[qi];
    Matrix3f Fprod = q->F[0];
    for (unsigned int jj = 1; jj < q->F.size(); jj++){
      Fprod = q->F[jj] * Fprod;
    }
    //fine element idx
    int feidx = q->ei[0];
    MaterialQuad* mat = (MaterialQuad* )(m[0]->m[m[0]->me[feidx]]);
    E += q->w[level] * mat->e[0]->getEnergy(Fprod);
  }
  return vol* E;
}

std::vector<Vector3f> ElementMeshHier::getForce(int level, int eIdx)
{
  Element * e = m[level]->e[eIdx];
  std::vector<Vector3f> fe(e->nV());

  return fe;
}

std::vector<Vector3f> ElementMeshHier::getForce(int level)
{
  ElementMesh * mesh = m[level];
  std::vector<Vector3f> f(mesh->x.size(), Vector3f(0,0,0));
  for (unsigned int ii = 0; ii < mesh->e.size(); ii++){
    std::vector<Vector3f> fe = getForce(level, ii);
    Element * e = mesh->e[ii];
    for (int jj = 0; jj < e->nV(); jj++){
      f[e->at(jj)] += fe[jj];
    }
  }
  return f;
}

ElementMeshHier::~ElementMeshHier()
{
}
