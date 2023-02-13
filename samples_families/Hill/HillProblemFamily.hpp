#pragma once

#include "IOptProblemFamily.hpp"
#include "HillProblem.hpp"

// Hill problem family
class THillProblemFamily : public IOptProblemFamily
{
public:
  THillProblemFamily() : IOptProblemFamily()
  {
    for (int i = 0; i < NUM_HILL_PROBLEMS; i++)
    {
      pOptProblems.push_back(new THillProblem(i));
    }
  }
};
