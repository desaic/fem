#ifndef ELEMENT_MESH_UTIL_HPP
#define ELEMENT_MESH_UTIL_HPP

#include "cfgDefs.h"
class ElementMesh;

std::vector<int> topVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> botVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> leftVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> rightVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> frontVerts(ElementMesh * em, const std::vector<int> & gridSize);
std::vector<int> backVerts(ElementMesh * em, const std::vector<int> & gridSize);
  ///@brief eight corners of a grid for coarsening.
std::vector<int> cornerVerts(ElementMesh * em, const std::vector<int> & gridSize);

///@brief assuming element size 1.
Vector3d shapeFunGrad(int ii, const Vector3d & xx);

///@param size. size of a cube element.
Eigen::MatrixXd BMatrix(const Vector3d & xx, const Eigen::Vector3d & size);

Eigen::VectorXd hexStrain(const Eigen::VectorXd & x, const Eigen::VectorXd & X,
  const Eigen::Vector3d & xi);

#endif