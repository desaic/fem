#ifndef HOMOGENIZE_3D_HPP
#define HOMOGENIZE_3D_HPP

#include <Eigen/Dense>
#include <vector>

class ElementMesh;
struct PardisoState;
class Homogenize3D{

public:
  
  Homogenize3D() :em(0), dim(3), gridSize(3,0){}
  
  void init();
  void solve();

  int dim;
  ElementMesh * em;
  std::vector<int> gridSize;
  
  std::vector<int> d1, d2;
  ///@brief edge dof followed by face dof.
  ///x y z axis then x y z faces.
  std::vector<int> d3, d4;
  std::vector<int> m_I, m_J;

  ///@brief corner displacements
  Eigen::MatrixXf ufixed;
  ///@brief edge displacements followed by face displacements
  Eigen::MatrixXf wfixed;
  ///@brief map from dof in d3 to d4
  ///d4 = M43d3 + wfixed
  Eigen::SparseMatrix<float> M43;
  ///@brief values in stiffness matrix.
  std::vector<double> m_val;
  Eigen::MatrixXd G;
  std::vector< std::vector<double> > u;

  PardisoState * pardisoState;
};

#endif