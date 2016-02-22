#ifndef exDefs_h
#define exDefs_h

#define cfgScalar float

#include "Eigen/Sparse"
typedef Eigen::Matrix<cfgScalar, 2, 1> Vector2S;
typedef Eigen::Matrix<cfgScalar, 2, 2> Matrix2S;
typedef Eigen::Matrix<cfgScalar, 3, 1> Vector3S;
typedef Eigen::Matrix<cfgScalar, 3, 3> Matrix3S;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef Eigen::Matrix<cfgScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXS;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixEXd;

#define Vec(d) Eigen::Matrix<cfgScalar, d, 1>

#include <vector>

#define SAFE_DELETE(x) { \
  if (x!=NULL) { \
    delete (x); \
    (x)=NULL; }\
  }

#endif 




