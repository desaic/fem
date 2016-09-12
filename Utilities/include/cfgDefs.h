#ifndef exDefs_h
#define exDefs_h

#ifdef USE_DOUBLE_FOR_SCALAR
  #define cfgScalar  double
#else
  #define cfgScalar  float
#endif

#include "Eigen/Sparse"
typedef Eigen::Matrix<cfgScalar, 2, 1> Vector2S;
typedef Eigen::Matrix<cfgScalar, 2, 2> Matrix2S;
typedef Eigen::Matrix<cfgScalar, 3, 1> Vector3S;
typedef Eigen::Matrix<cfgScalar, 3, 3> Matrix3S;
typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef Eigen::Matrix<cfgScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXS;

#define Vec(d) Eigen::Matrix<cfgScalar, d, 1>

#include <vector>

#define SAFE_DELETE(x) { \
  if (x!=NULL) { \
    delete (x); \
    (x)=NULL; }\
  }

#define SAFE_DELETE_VEC(p) { for (int i=0; i<p.size(); i++) {if(p[i] != NULL) {delete (p[i]); (p[i])=NULL; p.clear();}} }

#endif 




