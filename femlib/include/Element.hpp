#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <vector>

#include <Eigen/Dense>

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
  Eigen::Matrix3f defGrad(Eigen::Vector3f p,
    const std::vector<Eigen::Vector3f> & X, const std::vector<Eigen::Vector3f> & x) const;

  ///@param p natural coordinates.
  ///@return interpolation weight for each dof.
  virtual std::vector<float> shapeFun(const Eigen::Vector3f & p) const = 0;

  ///@param ii index of basis function.
  ///@param xx point in natural coordinates.
  ///@param X global array of rest positions.
  virtual Eigen::Vector3f shapeFunGrad(int ii, const Eigen::Vector3f & xx,
    const std::vector<Eigen::Vector3f> & X) const = 0;

  virtual std::vector<std::array<int,2> > getEdges();
  
  virtual float getVol(const std::vector<Eigen::Vector3f> & X);

  virtual Eigen::Vector3f getDisp(const Eigen::Vector3f & p, const std::vector<Eigen::Vector3f> & X,
    const std::vector<Eigen::Vector3f>x);

  const std::vector<int> & getNodeIndices()const{ return n; }
private:

  ///@brief nodal indices
  std::vector<int> n;

};
#endif