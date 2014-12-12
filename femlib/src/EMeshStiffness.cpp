#include "Element.hpp"
#include "ElementMesh.hpp"
#include "Material.hpp"
#include "MatrixXd.hpp"

MatrixXd ElementMesh::getStiffness(int eIdx)
{
  MatrixXd K = m[me[eIdx]]->getStiffness(e[eIdx],this);
  return K;
}

MatrixXd ElementMesh::getStiffness()
{
  int matSize = 3 * (int)x.size();
  MatrixXd Kglobal(matSize,matSize);
  Kglobal.fill(0);
  for(unsigned int ii = 0;ii<e.size();ii++){
    MatrixXd K = getStiffness(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      int vj = e[ii]->at(jj);
      for(int kk = 0; kk<e[ii]->nV(); kk++){
        int vk = e[ii]->at(kk);    
        addSubMat(K,Kglobal, 3*jj, 3*kk, 3*vj, 3*vk, 3, 3);
      }
    }
  }
  return Kglobal;
}


void getEleX(int ii, const ElementMesh * m, std::vector<Vector3f> &x)
{
  Element * ele = m->e[ii];
  x.resize(ele->nV());
  for (int jj = 0; jj<ele->nV(); jj++){
    x[jj] = m->x[ele->at(jj)];
  }
}

void setEleX(int ii, ElementMesh * m, const std::vector<Vector3f> &x)
{
  Element * ele = m->e[ii];
  for (int jj = 0; jj<ele->nV(); jj++){
    m->x[ele->at(jj)] = x[jj];
  }
}
