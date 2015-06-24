#ifndef ELEMENTHEX2D_HPP
#define ELEMENTHEX2D_HPP
#include "Element2D.h"

#include <Eigen\Dense>

class ElementHex2D:public Element2D
{
public:
  ElementHex2D();
  ElementHex2D(const ElementHex2D & e);
  int nV() const override{ return 4; }

  int nF() const override{ return 6; }

  std::vector<float> shapeFun(const Vector2f & p)const;

  Vector2f shapeFunGrad(int ii, const Vector2f & xx,
    const std::vector<Vector2f> & X) const;

  float getVol(const std::vector<Vector2f> & X);

  void setThickness(float iThickness) {m_thickness = iThickness;};
  float getThickness() {return m_thickness;};

private:
  float m_thickness;
};
#endif