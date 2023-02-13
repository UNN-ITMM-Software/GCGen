#pragma once

#include "grishagin_function.hpp"
#include "IOptProblemFamily.hpp"

class TGrishaginProblemFamily : public IOptProblemFamily
{
public:
  TGrishaginProblemFamily() : IOptProblemFamily()
  {
    for (int i = 1; i <= 100; i++)
      pOptProblems.push_back(new TGrishaginProblem(i));
  }
};