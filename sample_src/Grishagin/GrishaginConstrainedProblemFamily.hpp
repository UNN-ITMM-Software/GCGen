#pragma once

#include "GrishaginConstrainedProblem.hpp"
#include "IConstrainedOptProblemFamily.hpp"

class TGrishaginConstrainedProblemFamily : public IConstrainedOptProblemFamily
{
public:
  TGrishaginConstrainedProblemFamily(EConstrainedProblemType problemType = cptInFeasibleDomain,
    double fraction = 0.5, int activeConstrNum = 0) : IConstrainedOptProblemFamily()
  {
    for (int i = 1; i <= 100; i++)
      pOptProblems.push_back(new
        GrishaginConstrainedProblem(problemType, fraction, activeConstrNum,  i));
  }
};