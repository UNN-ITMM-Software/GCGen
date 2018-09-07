#pragma once

#include "GKLSProblem.hpp"
#include "IOptProblemFamily.hpp"

class TGKLSProblemFamily : public IOptProblemFamily
{
public:
  TGKLSProblemFamily(int dim = 2, GKLSClass type = Simple, GKLSFuncionType functionType = TD) : IOptProblemFamily()
  {
    for (int i = 1; i <= 100; i++)
      pOptProblems.push_back(new TGKLSProblem(i,dim, type, functionType));
  }
};