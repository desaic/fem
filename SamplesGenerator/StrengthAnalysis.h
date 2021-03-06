#ifndef StrengthAnalysis_h
#define StrengthAnalysis_h

#include "cfgDefs.h"

class ElementRegGrid2D;
class ElementRegGrid;
class MaterialQuad2D;
class MaterialQuad;
struct PardisoState;

class StrengthAnalysis
{
public:

  enum StructureType
  {
    General,
    Cubic,
    Orthotropic
  };

  StrengthAnalysis();
  ~StrengthAnalysis();

  bool run2D(int N[2], StructureType iType, std::vector<MaterialQuad2D> &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments,  const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensors, 
             int NbBlockRep, int NbSubdiv, std::vector<cfgScalar> &oStrengths);
  bool run3D(int N[3], StructureType iType, std::vector<MaterialQuad> &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensors, 
                            int iNbBlockRep, int iNbSubdiv, std::vector<cfgScalar> &oStrengths);

private:
  ElementRegGrid2D * createPhysicalSystem(int n[2], std::vector<MaterialQuad2D> &iMaterials, const std::vector<int> &iMatAssignments);
  ElementRegGrid * createPhysicalSystem(int n[3], std::vector<MaterialQuad> &iMaterials, const std::vector<int> &iMatAssignments);

  void initPardiso(ElementRegGrid2D *iPhysicalSystem);
  void initPardiso(ElementRegGrid *iPhysicalSystem);
  void getStiffnessSparse(ElementRegGrid2D * iPhysicalSystem, std::vector<double> &oValues, MatrixXS &oK);
  void getStiffnessSparse(ElementRegGrid* iPhysicalSystem, std::vector<double> &oValues, MatrixXS &oK);
  cfgScalar computeStrength(ElementRegGrid2D * iPhysicalSystem, MatrixXS &iElasticityTensor, const std::vector<cfgScalar>  &iYoungsModuliCoarse, const std::vector<cfgScalar> &iYoungsModulusBaseMaterials);
  cfgScalar computeStrength(ElementRegGrid* iPhysicalSystem, MatrixXS &iElasticityTensor, const std::vector<cfgScalar>  &iYoungsModuliCoarse, const std::vector<cfgScalar> &iYoungsModulusBaseMaterials);

private:
  ElementRegGrid2D * m_physicalSystem2D;
  ElementRegGrid * m_physicalSystem;
  std::vector<int> m_I, m_J;
  PardisoState * m_pardisoState; 
  std::vector<MaterialQuad2D> m_baseMaterials2D;
  std::vector<MaterialQuad> m_baseMaterials3D;
  int m_nbBlockRep;
  int m_nbSubdiv;
  MatrixXS m_K0[2];
  int m_dim;

  bool m_init;
  };

#endif





