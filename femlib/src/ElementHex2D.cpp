#include "ElementHex2D.h"

const int sw[4][2] =
{{-1,-1},
 {-1, 1}, 
 { 1,-1},
 { 1, 1},
};

///@brief face normals
const int facen[6][3]={
	{-1, 0, 0},
	{ 1, 0, 0},
	{ 0,-1, 0},
	{ 0, 1, 0},
	{ 0, 0,-1},
	{ 0, 0, 1}
};

ElementHex2D::ElementHex2D():Element2D(4)
{
  m_thickness = 1.f;
}

ElementHex2D::ElementHex2D(const ElementHex2D & e) :Element2D(e)
{
  m_thickness = 1.f;
}

ElementHex2D::ElementHex2D(const std::vector<int> &iNodes)
  :Element2D(iNodes)
{
  m_thickness = 1.f;
}

std::vector<cfgScalar>
ElementHex2D::shapeFun(const Vector2S & p)const
{
  std::vector<cfgScalar> weights(4);
  for(int ii = 0;ii<nV();ii++){
    weights[ii] = ((cfgScalar)1/(cfgScalar)4) * (1+sw[ii][0]*p[0]) *(1+sw[ii][1]*p[1]);
  }
  return weights;
}

Vector2S ElementHex2D::shapeFunGrad(int ii, const Vector2S & xx, const std::vector<Vector2S> & X) const
{
  Vector2S size=4*(X[at(3)] - X[at(0)]);
  Vector2S grad;
  size[0] = 1.0f/(size[0]);
  size[1] = 1.0f/(size[1]);

  grad[0] = sw[ii][0] * size[0] * (1 + sw[ii][1] * xx[1]);
  grad[1] = sw[ii][1] * size[1] * (1 + sw[ii][0] * xx[0]) ;

  return grad;
}

cfgScalar ElementHex2D::getVol(const std::vector<Vector2S> & X)
{
  Vector2S size = X[at(3)] - X[at(0)];
  return size[0] * size[1] * m_thickness;
}
