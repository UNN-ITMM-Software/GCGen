#pragma once

#include "IGeneralOptProblem.hpp"

/// Base class for the problem series
class IGeneralOptProblemFamily
{
protected:
  // pointers to the problems; problems must be generated in the class constructor
  // number of problems  = size of vector pOptProblems
  vector<IGeneralOptProblem*> pOptProblems;

  // The constructor is hidden so that one can not create an object of this class
  // It is required to declare a child, move the constructor to public and implement it
  IGeneralOptProblemFamily() {}
public:
  int GetFamilySize() const;
  ~IGeneralOptProblemFamily();
};
