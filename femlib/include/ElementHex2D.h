#ifndef ELEMENTHEX2D_HPP
#define ELEMENTHEX2D_HPP
#include "Element2D.h"

#include <Eigen/Dense>

class ElementHex2D:public Element2D
{
public:
  ElementHex2D();
  ElementHex2D(const ElementHex2D & e);
  ElementHex2D(const std::vector<int> &iNodes);

  int nV() const override{ return 4; }

  int nF() const override{ return 6; }

  std::vector<cfgScalar> shapeFun(const Vector2S & p)const;

  Vector2S shapeFunGrad(int ii, const Vector2S & xx,
    const std::vector<Vector2S> & X) const;

  cfgScalar getVol(const std::vector<Vector2S> & X);

  void setThickness(cfgScalar iThickness) {m_thickness = iThickness;};
  cfgScalar getThickness() {return m_thickness;};

private:
  cfgScalar m_thickness;
};
#endif