#include "FamilyGenerator.h"

#include "FamilyGeneratorImpl.h"

FamilyGenerator * FamilyGenerator::createOperator()
{
  FamilyGeneratorImpl * op = new FamilyGeneratorImpl();
  return op;
}
