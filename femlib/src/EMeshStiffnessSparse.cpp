#include "ElementMesh.hpp"
#include "Element.hpp"
#include "Eigen/Sparse"
#include "MatrixX.hpp"

typedef Eigen::Triplet<float> Tripletf;

void ElementMesh::getStiffnessSparse(std::vector<float> &val, bool trig, bool constrained)
{
  int N = 3* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXf K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            float val = K(3 * jj + dim1, 3 * kk + dim2);
            if (constrained){
              if (fixed[vk] || fixed[vj]){
                val = 0;
                if (vj == vk && dim1 == dim2){
                  val = 100;
                }
              }
            }
            //if (vk == vj && dim1 == dim2){
            //  val += 1;
            //}
            Tripletf triple(3 * vj + dim1, 3 * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     val.push_back(it.value());
   }
  }
}

void ElementMesh::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,
   bool trig)
{
  int N = 3* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);

  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    for(int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig&&(3*vk+dim2 > 3*vj+dim1)){
              continue;
            }
            Tripletf triple(3*vj+dim1,3*vk+dim2,1);
            coef.push_back(triple);
          }
        }
      }
    }
  }

  Ksparse.setFromTriplets(coef.begin(), coef.end());

  I.push_back(0);
  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     J.push_back(it.row());
   }
    I.push_back((int)J.size());
  }
}
