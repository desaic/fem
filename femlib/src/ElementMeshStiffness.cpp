
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
  MatrixXd Kglobal;
  Kglobal.allocate(matSize, matSize);
  for(unsigned int ii = 0;ii<e.size();ii++){
    MatrixXd K = getStiffness(ii);
    int nDof = 3* e[ii]->nV();
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      int vj = e[ii]->at(jj);
      for(int kk = 0; kk<e[ii]->nV(); kk++){
        int vk = e[ii]->at(kk);    
        copySubMat(K,Kglobal, 3*jj, 3*kk, 3*vj, 3*vk, 3, 3);
      }
    }
  }
  return Kglobal;
}

