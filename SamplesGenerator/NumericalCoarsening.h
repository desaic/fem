#ifndef NumericalCoarsening_h
#define NumericalCoarsening_h

#include "cfgDefs.h"
class ElementRegGrid2D;
class ElementRegGrid;
class MaterialQuad2D;
class MaterialQuad;
class Stepper2D;
class Stepper;
class Vector3f;
class PardisoState;

class NumericalCoarsening
{
public:
  enum StructureType
  {
    General,
    Cubic,
    Orthotropic
  };

  NumericalCoarsening();
  ~NumericalCoarsening();

  void init(int N[3], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad> &iBaseMaterials);

  void computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials, 
                                                     StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);
  void computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials, 
                                                     StructureType iType, std::vector<cfgScalar> &iBaseMaterialDensities, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);

  void computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[3], StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);
  void computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[3], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad> &iBaseMaterials, 
                                                     StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);

  MatrixXS computeCoarsenedElasticityTensor(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials);
  MatrixXS computeCoarsenedElasticityTensor(const std::vector<int> &iMaterials, int N[3], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad> &iBaseMaterials);

private:
  ElementRegGrid2D * createPhysicalSystem(int n[2], std::vector<MaterialQuad2D> &iMaterials, const std::vector<int> &iMatAssignments);
  ElementRegGrid * createPhysicalSystem(int n[3], std::vector<MaterialQuad> &iMaterials, const std::vector<int> &iMatAssignments);

  Stepper2D * createStepper(ElementRegGrid2D * physicalSystem);
  Stepper* createStepper(ElementRegGrid * physicalSystem);
  void addExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, Vector2S &iForceMagnitude);
  void addExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, Vector3S &iForceMagnitude);
  void setExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, ElementRegGrid2D *iPhysicalSystem);
  void setExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, cfgScalar iFz, ElementRegGrid *iPhysicalSystem);
  void getExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, Vector2S &iForceMagnitude, std::vector<cfgScalar> &oForces);
  void getExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, Vector3S &iForceMagnitude, std::vector<cfgScalar> &oForces);
  void getExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, ElementRegGrid2D *iPhysicalSystem, std::vector<float> &oForces);
  void getExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, cfgScalar iFz, ElementRegGrid *iPhysicalSystem, std::vector<float> &oForces);

  void computeHarmonicDisplacements(ElementRegGrid2D * iPhysicalSystem, Stepper2D * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements);
  void computeHarmonicDisplacements(ElementRegGrid * iPhysicalSystem, Stepper * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements);

  void computeHarmonicDisplacements(ElementRegGrid* iPhysicalSystem, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements);

  void getStiffnessSparse(ElementRegGrid * iPhysicalSystem, std::vector<double> &oValues);

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

  bool m_init;
};

#endif





