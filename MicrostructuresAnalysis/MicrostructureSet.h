#ifndef MicrostructureSet_h
#define MicrostructureSet_h

#include "cfgDefs.h"

class MicrostructureSet
{
public:
  MicrostructureSet();
  ~MicrostructureSet();

  void setMicrostructures(const std::vector<std::vector<int> > *iMatAssignments, int iMicrostructureSize[3], const std::vector<cfgScalar> *iMatParameters, int iParameterDim);

  int getParameterDim() const {return m_parameterDim;}
  int getNumStructures() const; 
  void getMicrostructureSize(int &oNx, int &oNy, int &oNz) const {oNx=m_microstructureSize[0]; oNy=m_microstructureSize[1]; oNz=m_microstructureSize[2];} 

  const std::vector<cfgScalar> & getAllMatParameters() const;
  const std::vector<int> & getMatAssignment(int iIndex) const;

private:
  const std::vector<std::vector<int> > *m_matAssignments;
  int m_microstructureSize[3];

  const std::vector<cfgScalar> *m_matParameters;
  int m_parameterDim;

  std::vector<int> m_voidVectorInt;
  std::vector<cfgScalar> m_voidVectorScalar;
};

#endif