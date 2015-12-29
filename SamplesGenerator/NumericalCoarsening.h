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

  void computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials, 
                                                     StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);
  void computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials, 
                                                     StructureType iType, std::vector<cfgScalar> &iBaseMaterialDensities, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters);

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

  void computeHarmonicDisplacements(ElementRegGrid2D * iPhysicalSystem, Stepper2D * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements);
  void computeHarmonicDisplacements(ElementRegGrid * iPhysicalSystem, Stepper * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements);
};

#endif





