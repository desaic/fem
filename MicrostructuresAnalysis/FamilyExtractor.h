#ifndef FamilyExtractor_h
#define FamilyExtractor_h

#include "cfgDefs.h"

class MicrostructureSet;

class FamilyExtractor
{
public:
  static FamilyExtractor * createOperator();
  virtual ~FamilyExtractor(){};

  virtual void setOption(int iOption){};
  virtual void setMicrostructures(MicrostructureSet *iMicrostructures, const std::vector<int> *iIndices=NULL) {};

  virtual bool step() = 0;
  virtual bool run() = 0;

  virtual const std::vector<int> & getMicrostructureIndices() = 0;
  virtual const std::vector<cfgScalar> & getReducedCoordinates() = 0;
  virtual int getReducedDimension() = 0;
};

#endif