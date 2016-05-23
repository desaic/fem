#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <vector>

#include <Eigen/Dense>
#include "cfgDefs.h"

class Element
{
public:
  
  Element(int _n=0);
  Element(const Element & e);
  Element(const std::vector<int> &inodes);

  ///@brief number of vertices
  virtual int nV() const;

  //@brief number of faces. Default returns 0 (not meaningful).
  virtual int nF() const;


  int & operator[](int idx);
  int operator[](int idx)const;
  int at(int idx)const;
  void resize(int size);

  ///@param X Reference vertex positions.
  ///@param x global array of Deformed vertex positions.
  Matrix3S defGrad(Vector3S p,
    const std::vector<Vector3S> & X, const std::vector<Vector3S> & x) const;

  ///@param p natural coordinates.
  ///@return interpolation weight for each dof.
  virtual std::vector<cfgScalar> shapeFun(const Vector3S & p) const = 0;

  ///@param ii index of basis function.
  ///@param xx point in natural coordinates.
  ///@param X global array of rest positions.
  virtual Vector3S shapeFunGrad(int ii, const Vector3S & xx,
    const std::vector<Vector3S> & X) const = 0;

  virtual std::vector<std::array<int,2> > getEdges();
  
  virtual cfgScalar getVol(const std::vector<Vector3S> & X);

  virtual Vector3S getDisp(const Vector3S & p, const std::vector<Vector3S> & X,
    const std::vector<Vector3S>x);

  const std::vector<int> & getNodeIndices()const{ return n; }

  //for rendering
  Vector3S color;
private:

  ///@brief nodal indices
  std::vector<int> n;

};
#endif