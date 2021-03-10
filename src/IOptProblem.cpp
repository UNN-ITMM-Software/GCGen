#include "IOptProblem.hpp"

// ------------------------------------------------------------------------------------------------
void IOptProblem::SetLipschitzConstant(double lipConst)
{
  IGeneralOptProblem::SetLipschitzConstant(0, lipConst);
}

// ------------------------------------------------------------------------------------------------
void IOptProblem::SetFunctionMax(vector<double> maxPoint, double maxValue)
{
  IGeneralOptProblem::SetFunctionMax(0, maxPoint, maxValue);
}

// ------------------------------------------------------------------------------------------------
IOptProblem::IOptProblem() : IGeneralOptProblem()
{
  mFunctionIndex = 0;
}

// ------------------------------------------------------------------------------------------------
IOptProblem::IOptProblem(int dim, vector<double> loBound, vector<double> upBound,
  vector<double> optPoint, double optVal, int probIndex)
  : IGeneralOptProblem(dim, loBound, upBound, probIndex)
{
  mOptimumPoint = optPoint;
  mOptimumValue = optVal;
  mFunctions.push_back({ mDimension, mOptimumPoint, mOptimumValue,{}, 0, -1, true, false,
    false, false });
  mFunctionNumber = mFunctions.size();
  mFunctionIndex = 0;
}

// ------------------------------------------------------------------------------------------------
vector<double> IOptProblem::GetOptimumPoint() const
{
  return IGeneralOptProblem::GetOptimumPoint();
}

// ------------------------------------------------------------------------------------------------
double IOptProblem::GetOptimumValue() const
{
  return IGeneralOptProblem::GetOptimumValue();
}

// ------------------------------------------------------------------------------------------------
bool IOptProblem::GetStatus(enum EOptFunctionParameter param) const
{
  return IGeneralOptProblem::GetStatus(mFunctionIndex, param);
}

// ------------------------------------------------------------------------------------------------
vector<double> IOptProblem::GetMaxPoint() const
{
  return IGeneralOptProblem::GetMaxPoint(mFunctionIndex);
}

// ------------------------------------------------------------------------------------------------
double IOptProblem::GetMaxValue() const
{
  return IGeneralOptProblem::GetMaxValue(mFunctionIndex);
}

// ------------------------------------------------------------------------------------------------
double IOptProblem::GetLipschitzConstant() const
{
  return IGeneralOptProblem::GetLipschitzConstant(mFunctionIndex);
}

// ------------------------------------------------------------------------------------------------
double IOptProblem::ComputeFunction(const vector<double>& y) const
{
  return Compute(mFunctionIndex, y);
}

// ------------------------------------------------------------------------------------------------
vector<double> IOptProblem::ComputeFunctionDerivatives(const vector<double>& y) const
{
  return ComputeDerivatives(mFunctionIndex, y);
}