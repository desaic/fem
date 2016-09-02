#include "ElementMesh.hpp"
#include "ElementHex.hpp"
#include "Material.hpp"
#include "Element.hpp"
#include "femError.hpp"
#include <iostream>
using namespace Eigen;
///@TODO only loads hex mesh
void ElementMesh::load(std::istream & in, cfgScalar scale)
{
  int nEle , nVert;
  std::string token;
  in>>token;
  in>>nVert;
  in>>token;
  in>>nEle;
  X.resize(nVert);
  e.resize(nEle);
  std::cout<<nVert<<" vertices\n"<<nEle<<" elements.\n";
  Vector3S minPos(0,0,0);
  for(int ii =0 ; ii<nVert; ii++){
    in>>X[ii][0]>>X[ii][1]>>X[ii][2];
    X[ii] = scale * X[ii];
    for(int jj =0 ; jj<3; jj++){
      if(minPos[jj]>X[ii][jj]){
        minPos[jj] = X[ii][jj];
      }
    }
  }
  std::cout<<minPos[0]<<" "<<minPos[1]<<" "<<minPos[2]<<"\n";
  for(int ii =0 ; ii<nVert; ii++){
    X[ii] -= minPos;
//    X[ii][1] += 0.001;
  }

  for(int ii =0; ii<nEle; ii++){
    int nV;
    in>>nV;
    Element * ele = 0;
    if(nV==8){
      ele = new ElementHex();
    }
    e[ii] = ele;
    for(int jj =0 ; jj<nV; jj++){
      in>>(*ele)[jj];
    }
  }
  e.resize(nEle);
  initArrays();
}


void ElementMesh::initArrays()
{
  x=X;
  v.resize(X.size(), Vector3S(0,0,0));
  fe.resize(X.size(), Vector3S(0,0,0));
  fixed.resize(X.size(),false);
  me.resize(e.size(),0);
}

void ElementMesh::addMaterial(Material*_m)
{
  m.push_back(_m);
  _m->init(this);
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

cfgScalar ElementMesh::getEnergy()
{
  cfgScalar ene = 0;
  for(unsigned int ii = 0;ii<e.size();ii++){
    ene += getEnergy(ii);
    if(fem_error){
      return -1;
    }
  }
  //energy from external forces
  for(unsigned int ii = 0;ii<fe.size();ii++){
    ene -= fe[ii].dot( x[ii]);
  }
  return ene;
}

cfgScalar ElementMesh::getEnergy(int eIdx)
{
  return m[me[eIdx]]->getEnergy(e[eIdx], this);
}

std::vector<Vector3S> ElementMesh::getForce(int eIdx)
{
  return m[me[eIdx]]->getForce(e[eIdx], this);
}

std::vector<Vector3S> ElementMesh::getForce()
{
  std::vector<Vector3S> force(x.size(), Vector3S::Zero());
  for(unsigned int ii = 0;ii<e.size();ii++){
    std::vector<Vector3S> fele = getForce(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      force[ e[ii]->at(jj) ] += fele[jj];
    }
  }
  for(unsigned int ii= 0;ii<fe.size();ii++){
    force[ii] += fe[ii];
//    std::cout<<fe[ii][1]<<"\n";
  }
  return force;
}

void ElementMesh::computeMass()
{
  mass.resize(x.size(), 0);
  for(unsigned int ii =0 ; ii<e.size(); ii++){
    cfgScalar size = X[e[ii]->at(7)][0] - X[e[ii]->at(0)][0];
    cfgScalar vol = size*size*size;
    cfgScalar nodeMass = 0.125 * vol * density;
    for(int jj =0 ; jj<e[ii]->nV(); jj++){
      mass[e[ii]->at(jj)] += nodeMass;
//      std::cout<<"m: "<<mass[e[ii]->at(jj)]<<"\n";
    }
  }
}

cfgScalar ElementMesh::eleSize()
{
  return X[e[0]->at(7)][0] - X[e[0]->at(0)][0];
}

ElementMesh::ElementMesh():dt(0.01),G(Vector3S(0,-9.8,0)), density(1000),
forceDrawingScale(1)
{}

ElementMesh::~ElementMesh()
{
  for(unsigned int ii = 0; ii<e.size();ii++){
    delete e[ii];
  }
}
