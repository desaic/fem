#include "ElementMesh.hpp"
#include "Material.hpp"
#include "Element.hpp"
#include "femError.hpp"
#include <iostream>

void ElementMesh::initArrays()
{
  x=X;
  fe.resize(X.size());
  fixed.resize(X.size());
  me.resize(e.size());
}

int ElementMesh::check()
{
  if(fe.size()!=x.size()){
    std::cout<<"external force and dof size differ\n";
    return -1;
  }
  if(e.size()!=me.size()){
    std::cout<<"material assignment and num elements differ\n";
    return -1;
  }
  for(unsigned int ii = 0;ii<me.size();ii++){
    if(me[ii]<0 || me[ii]>=m.size()){
      std::cout<<"Material index out of bounds\n";
      std::cout<<"ele: "<<ii<<", mat: "<<me[ii]<<"\n";
      return -1;
    }
  }
    
  if(fixed.size() != x.size()){
    std::cout<<"fixed array differ in size to vertex array \n";
    return -1;
  }
  return 0;
}

float ElementMesh::getEnergy()
{
  float ene = 0;
  for(unsigned int ii = 0;ii<e.size();ii++){
    ene += getEnergy(ii);
    if(fem_error){
      return -1;
    }
  }
  //energy from external forces
  for(unsigned int ii = 0;ii<fe.size();ii++){
    ene -= Vector3f::dot(fe[ii], x[ii]);
  }
  return ene;
}

float ElementMesh::getEnergy(int eIdx)
{
  return m[me[eIdx]]->getEnergy(e[eIdx], this);
}

std::vector<Vector3f> ElementMesh::getForce(int eIdx)
{
  return m[me[eIdx]]->getForce(e[eIdx], this);
}

std::vector<Vector3f> ElementMesh::getForce()
{
  std::vector<Vector3f> force(x.size());
  for(unsigned int ii = 0;ii<e.size();ii++){
    std::vector<Vector3f> fele = getForce(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      force[ e[ii]->at(jj) ] += fele[jj];
    }
  }
  for(unsigned int ii= 0;ii<fe.size();ii++){
    force[ii] += fe[ii];
  }
  return force;
}

ElementMesh::ElementMesh():u(0)
{}

ElementMesh::~ElementMesh()
{
  for(unsigned int ii = 0; ii<e.size();ii++){
    delete e[ii];
  }
}
