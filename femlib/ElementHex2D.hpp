#ifndef ELEMENTHEX2D_HPP
#define ELEMENTHEX2D_HPP
#include "Element2D.hpp"

#include <Eigen\Dense>

class ElementHex2D:public Element2D
{
public:
  ElementHex2D();
  ElementHex2D(const ElementHex2D & e);
  int nV() const override{ return 8; }

  int nF() const override{ return 6; }


  ///@brief natural coordinate for a point in reference space
  Vector2f natCoord(const Vector2f & p, const std::vector<Vector2f> & X);
  
  std::vector<float> shapeFun(const Vector2f & p)const;

  Vector2f shapeFunGrad(int ii, const Vector2f & xx,
    const std::vector<Vector2f> & X) const;

  std::vector<std::array<int,2> > getEdges();
  float getVol(const std::vector<Vector2f> & X);

  Vector2f facePt(int fi, const Vector2f & x);
  ///@param fi face index.
  Eigen::MatrixXf NMatrix(int fi);
  ///@brief shape function matrix.
  Eigen::MatrixXf HMatrix(const Vector2f & xx);
  //specific for linear strain (F+F^T).
  Eigen::MatrixXf BMatrix(const Vector2f & xx, const std::vector<Vector2f>X);
};
#endif