#ifndef Tensor_h
#define Tensor_h

#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"

class Tensor
{
public:
  // A_ijkl = iA[k][l][i][j]; B_abkl = iB[k][l][a][b]; C_ijkl = A_ijab B_abkl
  static void doubleContraction(const std::vector<std::vector<Matrix2S> > &iA, const std::vector<std::vector<Matrix2S> > &iB, std::vector<std::vector<Matrix2S> > &oC);
  static void add(const std::vector<std::vector<Matrix2S> > &iA, std::vector<std::vector<Matrix2S> > &ioB);

  static MatrixXS toMatrixRepresentation(const std::vector<std::vector<Matrix2S> > &iA);
  static MatrixXS toVoigtRepresentation(const MatrixXS &iA);
  static MatrixXS toVoigtRepresentation(const std::vector<std::vector<Matrix2S> > &iA);   // A_ijkl = iA[k][l](i,j)
};

#endif