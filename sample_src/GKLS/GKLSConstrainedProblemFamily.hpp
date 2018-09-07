#pragma once

#include "GKLSConstrainedProblem.hpp"
#include "IConstrainedOptProblemFamily.hpp"

class TGKLSConstrainedProblemFamily : public IConstrainedOptProblemFamily
{
public:
  TGKLSConstrainedProblemFamily(EConstrainedProblemType problemType = cptInFeasibleDomain,
    double fraction = 0.5, int activeConstrNum = 0, int dim = 2,
    GKLSClass type = Simple, GKLSFuncionType functionType = TD) : IConstrainedOptProblemFamily()
  {
    for (int i = 1; i <= 100; i++)
      pOptProblems.push_back(new
        TGKLSConstrainedProblem(problemType, fraction, activeConstrNum,
          i, dim, type, functionType));
  }
};