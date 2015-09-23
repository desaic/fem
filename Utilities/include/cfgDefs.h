#ifndef exDefs_h
#define exDefs_h

#define cfgScalar float

#include "Eigen/Sparse"
typedef Eigen::Matrix<cfgScalar, 2, 1> Vector2S;
typedef Eigen::Matrix<cfgScalar, 2, 2> Matrix2S;
typedef Eigen::Matrix<cfgScalar, 3, 1> Vector3S;
typedef Eigen::Matrix<cfgScalar, 3, 3> Matrix3S;
typedef Eigen::Matrix<cfgScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXS;

#include <vector>

#endif 

