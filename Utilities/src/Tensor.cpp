#include "Tensor.h"

#include <assert.h>
#include <iostream>

// A_ijkl = iA[k][l](i,j); B_klab = iB[a][b](i,j); C_ijab = A_ijkl B_klab
void Tensor::doubleContraction(const std::vector<std::vector<Matrix2S> > &iA, const std::vector<std::vector<Matrix2S> > &iB, std::vector<std::vector<Matrix2S> > &oC)
{
  oC.clear();

  int n_a = (int)iB.size();
  oC.resize(n_a);
  for (int a=0; a<n_a; a++)
  {
    int n_b=(int)iB[a].size();
    for (int b=0; b<n_b; b++)
    {
      Matrix2S Cijab;
      int n_i=2;
      for (int i=0; i<n_i; i++)
      {
        int n_j=2;
        for (int j=0; j<n_j; j++)
        {
          Cijab(i,j) = 0;
          int n_k=2;
          for (int k=0; k<n_k; k++)
          {
            int n_l=2;
            for (int l=0; l<n_l; l++)
            {
              Cijab(i,j) += iA[k][l](i,j);
            }
          }
        }
      }
      oC[a].push_back(Cijab);
    }
  }
}

void Tensor::add(const std::vector<std::vector<Matrix2S> > &iA, std::vector<std::vector<Matrix2S> > &ioB)
{
  assert(iA.size()==ioB.size());
  int k, n_k = (int)iA.size();
  for (k=0; k<n_k; k++)
  {
    int l, n_l = (int)iA[k].size();
    assert(ioB[k].size()==n_l);
    for (l=0; l<n_l; l++)
    {
      ioB[k][l] += iA[k][l];
    }
  }
}

// A_ijkl = iA[k][l](i,j)
MatrixXS Tensor::toMatrixRepresentation(const std::vector<std::vector<Matrix2S> > &iA)
{
  int dim=2;
  int n_k = iA.size();
  int n_l = iA[0].size();

  MatrixXS M = MatrixXS::Zero(2*n_k,2*n_l);
  for (int k=0; k<n_k; k++)
  {
    for (int l=0; l<n_l; l++)
    {
      for (int i=0; i<2; i++)
      {
        for (int j=0; j<2; j++)
        {
          cfgScalar  A_ijkl = iA[k][l](i,j);
          M(2*i+j, 2*k+l) = A_ijkl;
        }
      }
    }
  }
  return M;
}

// A_ijkl = iA(k,l)[i][j]
MatrixXS Tensor::toVoigtRepresentation(const MatrixXS &iA)
{
  MatrixXS M = MatrixXS::Zero(3,3);
  int ind[3] = {0, 2, 3};
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      M(i, j) = iA(ind[i], ind[j]);
    }
  }
  return M;
}

// A_ijkl = iA[k][l](i,j)
MatrixXS Tensor::toVoigtRepresentation(const std::vector<std::vector<Matrix2S> > &iA)
{
  MatrixXS M1 = Tensor::toMatrixRepresentation(iA);
  std::cout << M1 << std::endl;


  MatrixXS M2 = Tensor::toVoigtRepresentation(M1);
  return M2;
}


