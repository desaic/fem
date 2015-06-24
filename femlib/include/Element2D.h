#ifndef ELEMENT2D_HPP
#define ELEMENT2D_HPP

#include <array>
#include <vector>

#include "vecmath.h"

class Element2D
{
public:
  
  Element2D(int _n=0);
  Element2D(const Element2D & e);
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
  Matrix2f defGrad(Vector2f p,
    const std::vector<Vector2f> & X, const std::vector<Vector2f> & x) const;

  ///@param p natural coordinates.
  ///@return interpolation weight for each dof.
  virtual std::vector<float> shapeFun(const Vector2f & p) const=0;

  ///@param ii index of basis function.
  ///@param xx point in natural coordinates.
  ///@param X global array of rest positions.
  virtual Vector2f shapeFunGrad(int ii, const Vector2f & xx,
    const std::vector<Vector2f> & X) const=0;

  virtual std::vector<std::array<int,2> > getEdges();
  
  virtual float getVol(const std::vector<Vector2f> & X);

  virtual Vector2f getDisp(const Vector2f & p, const std::vector<Vector2f> & X,
    const std::vector<Vector2f>x);

  const std::vector<int> & getNodeIndices()const{ return n; }
private:

  ///@brief nodal indices
  std::vector<int> n;

};
#endif