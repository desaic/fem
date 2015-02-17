#include "ElementHier.hpp"
#include "ElementMeshHier.hpp"
#include "ElementMesh.hpp"
#include "ElementRegGrid.hpp"
#include "Element.hpp"
#include "MaterialQuad.hpp"
#include <iostream>
const float eps = 1e-6;
///@brief construct nLevel additional levels given m[0].
void ElementMeshHier::buildHier(int nLevel)
{
  if (m.size() == 0){
    std::cout << "Error ElementMeshHier coarsen. Need fine mesh.\n";
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
}

ElementMeshHier::ElementMeshHier()
{
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


  return E;
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
