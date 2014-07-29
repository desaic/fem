#include "ElementMesh.hpp"
#include "Material.hpp"
#include <iostream>
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

  if(bd.size()!=x.size()){
    std::cout<<"number of bounds and num dof differ\n";
    return -1;
  }
  for(unsigned int ii = 0;ii<bd.size();ii++){
    for(unsigned int jj =0;jj<3;jj++){
      if(bd[ii].u[jj]<x[ii][jj]
      ||bd[ii].l[jj]>x[ii][jj]){
        std::cout<<"Initial configurature violates bounds\n";
        std::cout<<"vert: "<<ii<<", dim:"<<jj<<"\n";
        return -1;
      }
    }
  }
  return 0;
}

float ElementMesh::getEnergy()
{
  float ene = 0;
  for(unsigned int ii = 0;ii<e.size();ii++){
    ene += getEnergy(ii);
  }
  //energy from external forces
  for(unsigned int ii = 0;ii<fe.size();ii++){
    ene += Vector3f::dot(fe[ii], x[ii]);
  }
  return ene;
}

float ElementMesh::getEnergy(int eIdx)
{
  return m[me[eIdx]]->getEnergy(e[eIdx], this);
}

ElementMesh::ElementMesh()
{}