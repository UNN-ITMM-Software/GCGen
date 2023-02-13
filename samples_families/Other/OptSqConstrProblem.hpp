#pragma once

#include "IConstrainedOptProblem.hpp"

// Example of class for constrained problem
class TOptSqConstrProblem : public IConstrainedOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 0;
    switch (index)
    {
    case 0:
      for (int i = 0; i < mDimension; i++)
        res += y[i] * y[i];
      break;
    case 1:
      for (int i = 0; i < mDimension; i++)
        res += (y[i] - 0.5) * (y[i] - 0.5);
      res -= 4;
      break;
    default:
      throw string("Function index is out of range");
    }
    return res;
  }
public:
  TOptSqConstrProblem(int dim, vector<double> loBound, vector<double> upBound,
    vector<double> optPoint, double optVal,
    EConstrainedProblemType problemType = cptNormal, double fraction = 1, int activeConstrNum = 0)
    : IConstrainedOptProblem(dim, loBound, upBound, -1, problemType,
      fraction, activeConstrNum)
  {
    mOptimumPoint = optPoint;
    mOptimumValue = optVal;

    mFunctions[0].mMinimumPoint = optPoint;
    mFunctions[0].mMinimumValue = optVal;
    mFunctions[0].mIsMinimumKnown = true;

    vector<double> constrainOptimumPoint(mDimension);
    for (int i = 0; i < mDimension; i++)
      constrainOptimumPoint[i] = 0.5;
    mFunctions.push_back({ mDimension, constrainOptimumPoint, -4, {}, 0, -1,
      true, false, false, false });

    mFunctionNumber = mFunctions.size();

    mCriterionIndeces.push_back(0);

    mConstraintIndeces.push_back(1);

    InitProblem();
  }
};