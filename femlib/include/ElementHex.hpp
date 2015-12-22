#ifndef ELEMENTHEX_HPP
#define ELEMENTHEX_HPP
#include "Element.hpp"

#include <Eigen/Dense>
extern const int facen[][3];
extern const int faces[][4];
class ElementHex:public Element
{
public:
  ElementHex();
  ElementHex(const ElementHex & e);
  ElementHex(const std::vector<int> &iNodes);

  int nV() const override{ return 8; }

  int nF() const override{ return 6; }


  ///@brief natural coordinate for a point in reference space
  Eigen::Vector3f natCoord(const Eigen::Vector3f & p, const std::vector<Eigen::Vector3f> & X);
  
  std::vector<float> shapeFun(const Eigen::Vector3f & p)const;

  Eigen::Vector3f shapeFunGrad(int ii, const Eigen::Vector3f & xx,
    const std::vector<Eigen::Vector3f> & X) const;

  std::vector<std::array<int,2> > getEdges();
  float getVol(const std::vector<Eigen::Vector3f> & X);

  Eigen::Vector3f facePt(int fi, const Eigen::Vector3f & x);
  ///@param fi face index.
  Eigen::MatrixXf NMatrix(int fi);
  ///@brief shape function matrix.
  Eigen::MatrixXf HMatrix(const Eigen::Vector3f & xx);
  //specific for linear strain (F+F^T).
  Eigen::MatrixXf BMatrix(const Eigen::Vector3f & xx, const std::vector<Eigen::Vector3f>X);
};
#endif
