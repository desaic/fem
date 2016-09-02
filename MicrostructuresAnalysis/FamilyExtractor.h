#ifndef FamilyExtractor_h
#define FamilyExtractor_h

class MicrostructureSet;

class FamilyExtractor
{
public:
  static FamilyExtractor * createOperator();
  virtual ~FamilyExtractor(){};

  virtual void setMicrostructures(MicrostructureSet *iMicrostructures) {};

  virtual bool step() = 0;
  virtual bool run() = 0;
};

#endif