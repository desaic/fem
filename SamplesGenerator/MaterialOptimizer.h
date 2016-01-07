#ifndef MaterialOptimizer_h
#define MaterialOptimizer_h

#include "cfgDefs.h"

class MaterialQuad2D;
class MaterialQuad;
class RealFun;
class FEM2DFun;
class FEM3DFun;
class ElementRegGrid2D;

class MaterialOptimizer
{
public:
  enum StructureType
  {
    General,
    Cubic,
    Orthotropic
  };

  MaterialOptimizer();
  ~MaterialOptimizer();

  void setStructureType(StructureType iType);

  bool run2D(int N[2], const std::vector<MaterialQuad2D> &iBaseMaterials, const std::vector<int> &iMaterialAssignments, const std::vector<float> &iTargetParams, std::vector<std::vector<int> > &oNewMaterialAssignments);
  bool run3D(int N[3], const std::vector<MaterialQuad> &iBaseMaterials, const std::vector<int> &iMaterialAssignments, const std::vector<float> &iTargetParams, std::vector<std::vector<int> > &oNewMaterialAssignments);

private:

  ///@brief verify consistency of f and df using central differencing
  void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h);

  ///@brief finds a step size the decreases value of fun.
  ///@return <0 if cannot find a positive step size.
  ///@param x0 starting point.
  ///@param dir search direction.
  ///@param h appropriate step size if there is one. Always positive.
  int lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h);

  ///@param nSteps maximum number of steps.
  void gradientDescent(FEM2DFun * fun, Eigen::VectorXd & x0, int nSteps, std::vector<Eigen::VectorXd> &oParameters);
  void gradientDescent(FEM3DFun * fun, Eigen::VectorXd & x0, int nSteps, std::vector<Eigen::VectorXd> &oParameters);

  double infNorm(const Eigen::VectorXd & a);

  void computeTargetDisplacements(float iYoungModulus, float iPoissonsRatio, float iFx, float &odx, float &ody);
  void getExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, Vector2S &iForceMagnitude, std::vector<double> &oForces);

  Eigen::MatrixXd computeTargetStrains(int idim, float iYoungModulus, float iPoissonRatio, float iShearModulus);
  Eigen::MatrixXd computeTargetStrains(float iE_x, float iE_y, float iNu_xy, float iMu_xy);

private:
  StructureType m_structureType;

  };

#endif





