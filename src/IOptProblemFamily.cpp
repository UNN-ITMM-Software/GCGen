#include "IOptProblemFamily.hpp"

IOptProblem * IOptProblemFamily::operator[](int index)
{
  IOptProblem* p = static_cast<IOptProblem*>(pOptProblems[index]);
  return p;
}