#ifndef FamilyGenerator_h
#define FamilyGenerator_h

#include "cfgDefs.h"

class FamilyGenerator
{
public:
  static FamilyGenerator * createOperator();
  virtual ~FamilyGenerator(){};

  virtual void setMicrostructureSize(int iDim, int n[3]) = 0;

  virtual bool step() = 0;
  virtual bool run() = 0;

  virtual const std::vector<std::vector<int> > & getMicrostructures() = 0;
  virtual const std::vector<cfgScalar> & getParameters() = 0;
  virtual void getSize(int oN[3]) = 0;
};

#endif