#include "IConstrainedOptProblemFamily.hpp"

IConstrainedOptProblem* IConstrainedOptProblemFamily::operator[](int index)
{
  IConstrainedOptProblem* p = static_cast<IConstrainedOptProblem*>(pOptProblems[index]);
  return p;
}