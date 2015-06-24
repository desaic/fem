
#include "Element2D.h"
#include "ElementRegGrid2D.h"
#include "ElementHex2D.h"
#include "util.h"
#include <set>

//node index
#define IX(ii,jj) ( (ii) * NY + jj)

///vertex positions for single element
int oneEleV[4][2] =
{ {0,0},
  {0,1},
  {1,0},
  {1,1}
};

int ElementRegGrid2D::GetVertInd(int ii , int jj)
{
  int vidx =  ii * (ny+1) + jj;
  if(vertIdx.size()==0){
    return vidx;
  }
  return vertIdx[vidx];
}

int
ElementRegGrid2D::GetEleInd(int ii , int jj)const
{
  if(eleIdx.size()==0){
    return toArrInd(ii,jj);
  }
  return eleIdx[toArrInd(ii,jj)];
}

int
ElementRegGrid2D::GetEleInd(const Vector2f & p)
{
  Vector2f lp = p - origin;
  int ii = (int)( lp[0]/dx);
  int jj = (int)( lp[1]/dx);
//  ii = clamp(ii, 0, nx - 1);
//  jj = clamp(jj, 0, ny - 1);
  return GetEleInd(ii,jj);
}

int
ElementRegGrid2D::GetEleInd_clamp(const Vector2f & p)
{
  Vector3f lp = p - origin;
  int ii = (int)( lp[0]/dx);
  int jj = (int)( lp[1]/dx);
  ii = clamp(ii, 0, nx - 1);
  jj = clamp(jj, 0, ny - 1);
  return GetEleInd(ii,jj);
}

void
ElementRegGrid2D::rmEmptyVert()
{
  vertIdx.resize(X.size());
  std::fill(vertIdx.begin(), vertIdx.end(), -1);
  std::vector<Vector2f> newNodes;
  int vCnt =0;
  for(unsigned int ii = 0;ii<e.size(); ii++){
    Element2D * ele = e[ii];
    
    for(int jj = 0;jj<ele->nV();jj++){
      int vi = ele->at(jj);
      if(vertIdx[vi]<0){
        newNodes.push_back(X[vi]);
        vertIdx[vi]=vCnt;
        vCnt++;
      }
      vi = vertIdx[vi];
      (*ele)[jj] = vi;
    }
  }
  X=newNodes;
  x=X;
}

void
ElementRegGrid2D::deleteEle()
{
  for(auto ii = 0;ii< e.size();ii++){
    delete e[ii];
  }
}

void ElementRegGrid2D::allocate()
{
  for(auto ii = 0;ii< e.size();ii++){
    delete e[ii];
  }
  e.clear();
  X.clear();
  int maxn = std::max(nx, ny);
  dx = 1.0f / maxn;
  //number of vertices is 1 more than grid cells in each dimension
  int NY = ny + 1;
  for (int ii = 0; ii <= nx; ii++) {
    for (int jj = 0; jj <= ny; jj++) {
      X.push_back(dx * Vector2f((float)ii, (float)jj));
    }
  }
  for (int ii = 0; ii < nx; ii++) {
    for (int jj = 0; jj < ny; jj++) {
      std::vector<int> indices(4);
      ElementHex2D * ele = new ElementHex2D();
      for(int ll =0;ll<4;ll++){
        (*ele)[ll] = IX(ii+oneEleV[ll][0], jj+oneEleV[ll][1]);
      }
      e.push_back(ele);
    }
  }
  initArrays();
}

ElementRegGrid2D::ElementRegGrid2D(int _nx , int _ny )
  :nx(_nx),ny(_ny)
{
  if (nx>0 && ny>0){
    allocate();
  }
}

void
ElementRegGrid2D::resize(int _nx, int _ny)
{
  nx = _nx;
  ny = _ny;
  allocate();  
}

ElementRegGrid2D::~ElementRegGrid2D()
{
  //elements are deleted in base class ElementMesh.
}
