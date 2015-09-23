#ifndef ElementRegGrid2D2D_HPP_
#define ElementRegGrid2D2D_HPP_
#include "ElementMesh2D.h"
#include "vecmath.h"
class ElementRegGrid2D:public ElementMesh2D
{
public:
  ElementRegGrid2D(int _nx = 0 , int _ny = 0);
  void resize(int _nx, int _ny);
  void allocate();
  void deleteEle();

  ///@brief remove unused vertices
  void rmEmptyVert();
  
  int GetEleInd(const Vector2S & p);

  ///@brief 
  int GetEleInd(int ii , int jj)const;
  int GetEleInd_clamp(const Vector2S & p);
  int GetVertInd(int ii , int jj);
  int toArrInd(int ii , int jj)const{return ii * ny + jj;}

  virtual ~ElementRegGrid2D();
  std::vector<int> eleIdx;
  //map from  grid index to vertex index in array X and x.
  std::vector<int> vertIdx;
  int nx,ny;
  cfgScalar dx;

  Vector2S origin;
};
#endif /* ElementRegGrid2D_HPP_ */
