#include "FamilyExtractor.h"

#include "FamilyExtractorImpl.h"

FamilyExtractor * FamilyExtractor::createOperator()
{
  FamilyExtractorImpl * op = new FamilyExtractorImpl();
  return op;
}
