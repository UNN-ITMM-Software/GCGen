#include "IGeneralOptProblemFamily.hpp"

int IGeneralOptProblemFamily::GetFamilySize() const
{
  return pOptProblems.size();
}

IGeneralOptProblemFamily::~IGeneralOptProblemFamily()
{
  for (size_t i = 0; i < pOptProblems.size(); i++)
    delete pOptProblems[i];
}