#ifndef ELEMENT_COARSE_HPP
#define ELEMENT_COARSE_HPP

#include "ElementHex.hpp"
#include <vector>

class ElementRegGrid;

class ElementCoarse:public ElementHex
{
public:
  
  ElementCoarse();

  ///@brief Natural coordinates of 27 fine vertices
  static std::vector<Vector3f> Xfine;

  ///@brief index of fine element vertices in coarse element
  static std::vector<std::vector< int> > vifine;
  
  ///@brief index of fine elements in the fine mesh
  std::vector<int> fineEle;

  ///@param p natural coordinates
  Vector3f getDisp(Vector3f p);

  
};

///@bug grid assumed to have even dimensions for the moment.
std::vector<ElementCoarse *> coarsen(const ElementRegGrid & grid);

#endif
