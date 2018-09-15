#pragma once
#include "IConstrainedOptProblem.hpp"
#include "IGeneralOptProblemFamily.hpp"

class IConstrainedOptProblemFamily : public IGeneralOptProblemFamily
{
protected:
  // The constructor is hidden so that one can not create an object of this class
  // It is required to declare a child, move the constructor to public and implement it
  IConstrainedOptProblemFamily() {}
public:
  IConstrainedOptProblem* operator[](int index);
};
