#ifndef NumericalCoarseningUsingMultigrid_h
#define NumericalCoarseningUsingMultigrid_h

#include "cfgDefs.h"

class MaterialQuad;

#include "ElasticHexFEM.h"
#include "Homogenization.h"

class NumericalCoarseningUsingMultigrid
{
public:
  enum StructureType
  {
    General,
    Cubic,
    Orthotropic
  };

public:
  NumericalCoarseningUsingMultigrid();
  ~NumericalCoarseningUsingMultigrid();

  void setUseVariableDensity(bool iUseVariableDensity) {m_useVariableDensity = iUseVariableDensity;}

  void init(int N[3], int iNbBlockRep, int iNbSubdiv, const std::vector<MaterialQuad> &iBaseMaterials, StructureType iType);

  void computeCoarsenedElasticityTensorAndParameters3D(const std::vector<int> &iMaterials, int N[3], const std::vector<cfgScalar> &iBaseMaterialDensities, StructureType iType, 
                                                       std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);

private:

private:
  std::vector<MaterialQuad> m_baseMaterials3D;
  int m_nbBlockRep;
  int m_nbSubdiv;
  int m_dim;

  ElasticHexFEM<3> m_fem;
  Homogenization<3> * m_hom;
  std::vector<Eigen::VectorXd> m_u;
  std::vector<Eigen::MatrixXd> m_K0;
  bool m_useVariableDensity;
  cfgScalar m_densityMin;

  bool m_init;
};

#endif





