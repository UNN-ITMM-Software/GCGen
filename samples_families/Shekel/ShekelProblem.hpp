#pragma once

#include <vector>
#include <cmath>

#include "IOptProblem.hpp"

using std::vector;

#define NUM_SHEKEL_PROBLEMS 1000
#define NUM_SHEKEL_COEFF 10

extern double kShekel[NUM_SHEKEL_PROBLEMS][NUM_SHEKEL_COEFF];
extern double aShekel[NUM_SHEKEL_PROBLEMS][NUM_SHEKEL_COEFF];
extern double cShekel[NUM_SHEKEL_PROBLEMS][NUM_SHEKEL_COEFF];
extern double minShekel[NUM_SHEKEL_PROBLEMS][2];
extern double maxShekel[NUM_SHEKEL_PROBLEMS][2];
extern double lConstantShekel[NUM_SHEKEL_PROBLEMS];

// Shekel problem
class TShekelProblem : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 0.0;

    for (int i = 0; i < NUM_SHEKEL_COEFF; i++)
    {
      res = res - 1 / (kShekel[mProblemIndex][i] * pow(y[0] - aShekel[mProblemIndex][i], 2.0) +
        cShekel[mProblemIndex][i]);
    }
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 0.0;

    for (int i = 0; i < NUM_SHEKEL_COEFF; i++)
    {
      res = res + (2 * kShekel[mProblemIndex][i] * (y[0] - aShekel[mProblemIndex][i])) /
        pow((kShekel[mProblemIndex][i] * pow(y[0] - aShekel[mProblemIndex][i], 2.0) + cShekel[mFunctionIndex][i]), 2);
    }

    return{ res };
  }
public:
  TShekelProblem(int functionIndex = 0) : IOptProblem(1, { 0.0 }, { 10.0 },
  { minShekel[functionIndex][1] }, minShekel[functionIndex][0], functionIndex)
  {
    mProblemIndex = functionIndex;
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(lConstantShekel[functionIndex]);
    SetFunctionMax({ maxShekel[functionIndex][1] }, maxShekel[functionIndex][0]);
  }
};
