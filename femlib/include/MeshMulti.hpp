#ifndef MESHMULTI_HPP
#define MESHMULTI_HPP

#include "ElementRegGrid.hpp"

class ElementHex;

///@brief index of an element in a heirarchy.
struct HIndex
{
  ///@brief 0 is the fine level.
  int level;

  int idx;

};

class MeshMulti
{
public:

  std::vector<std::vector<ElementHex*> >eh;
  
  std::vector<std::vector<Vector3f> >xh;
  std::vector<std::vector<Vector3f> >Xh;

  ElementHex* operator[](const HIndex& ind){
    return eh[ind.level][ind.idx];
  }

  Vector3f getDisp(const HIndex & ind, const Vector3f & p);
};

#endif