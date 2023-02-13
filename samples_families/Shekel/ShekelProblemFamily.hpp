#pragma once

#include "IOptProblemFamily.hpp"
#include "ShekelProblem.hpp"

// Shekel problem family
class TShekelProblemFamily : public IOptProblemFamily
{
public:
  TShekelProblemFamily() : IOptProblemFamily()
  {
    for (int i = 0; i < NUM_SHEKEL_PROBLEMS; i++)
    {
      pOptProblems.push_back(new TShekelProblem(i));
    }
  }
};