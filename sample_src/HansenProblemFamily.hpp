#pragma once

#include "HansenProblem.hpp"
#include "IOptProblemFamily.hpp"


class THansenProblemFamily : public IOptProblemFamily
{
public:
  THansenProblemFamily() : IOptProblemFamily()
  {
    pOptProblems.push_back(new THansenProblem0());
    pOptProblems.push_back(new THansenProblem1());
    pOptProblems.push_back(new THansenProblem2());
    pOptProblems.push_back(new THansenProblem3());
    pOptProblems.push_back(new THansenProblem4());
    pOptProblems.push_back(new THansenProblem5());
    pOptProblems.push_back(new THansenProblem6());
    pOptProblems.push_back(new THansenProblem7());
    pOptProblems.push_back(new THansenProblem8());
    pOptProblems.push_back(new THansenProblem9());
    pOptProblems.push_back(new THansenProblem10());
    pOptProblems.push_back(new THansenProblem11());
    pOptProblems.push_back(new THansenProblem12());
    pOptProblems.push_back(new THansenProblem13());
    pOptProblems.push_back(new THansenProblem14());
    pOptProblems.push_back(new THansenProblem15());
    pOptProblems.push_back(new THansenProblem16());
    pOptProblems.push_back(new THansenProblem17());
    pOptProblems.push_back(new THansenProblem18());
    pOptProblems.push_back(new THansenProblem19());
  }
};