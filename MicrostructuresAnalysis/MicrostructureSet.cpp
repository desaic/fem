#include "MicrostructureSet.h"

#include <assert.h>

MicrostructureSet::MicrostructureSet()
{
  m_matAssignments = NULL;
  m_matParameters = NULL;

  m_microstructureSize[0] = 0;
  m_microstructureSize[1] = 0;
  m_microstructureSize[2] = 0;

  m_parameterDim = 0;
}

MicrostructureSet::~MicrostructureSet()
{
}

 void MicrostructureSet::setMicrostructures(const std::vector<std::vector<int>> *iMatAssignments, int iMicrostructureSize[3], const std::vector<cfgScalar> *iMatParameters, int iParameterDim)
 {
   m_matAssignments = iMatAssignments;
   for (int idim=0; idim<3; idim++)
   {
     m_microstructureSize[idim] = iMicrostructureSize[idim];
   }
   m_matParameters = iMatParameters;
   m_parameterDim = iParameterDim;
 }

 int MicrostructureSet::getNumStructures() const
 {
   if (m_matAssignments)
   {
     return (int)m_matAssignments->size();
   }
   else
   {
     return 0;
   }
 }

 const std::vector<cfgScalar> & MicrostructureSet::getAllMatParameters() const
 {
   if (m_matParameters)
   {
     return *m_matParameters;
   }
   else
   {
     return m_voidVectorScalar;
   }
 }

 const std::vector<int> & MicrostructureSet::getMatAssignment(int iIndex) const
 {
   if (m_matAssignments)
   {
     assert(iIndex>=0 && iIndex<m_matAssignments.size());
     return (*m_matAssignments)[iIndex];
   }
   else
   {
     m_voidVectorInt;
   }
 }

