#include "IGeneralOptProblem.hpp"

// ------------------------------------------------------------------------------------------------
bool IGeneralOptProblem::GetStatus(uint index, enum EOptFunctionParameter param) const
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  bool status = true;
  switch (param)
  {
  case ofpLipschitz:
    if (mFunctions[index].mIsLipschitzConstantKnown == false)
      status = false;
    break;
  case ofpMinimum:
    if (mFunctions[index].mIsMinimumKnown == false)
      status = false;
    break;
  case ofpMaximum:
    if (mFunctions[index].mIsMaximumKnown == false)
      status = false;
    break;
  case ofpDerivatives:
    if (mFunctions[index].mIsDerivativesKnown == false)
      status = false;
    break;
  default:
    status = false;
  }
  return status;
}

// ------------------------------------------------------------------------------------------------
void IGeneralOptProblem::SetLipschitzConstant(uint index, double lipConst)
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  mFunctions[index].mLipschitzConstant = lipConst;
  mFunctions[index].mIsLipschitzConstantKnown = true;
}

// ------------------------------------------------------------------------------------------------
double IGeneralOptProblem::GetLipschitzConstant(uint index) const
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  if (mFunctions[index].mIsLipschitzConstantKnown == false)
    throw string("Lipschitz constant is unknown");
  return mFunctions[index].mLipschitzConstant;
}

// ------------------------------------------------------------------------------------------------
void IGeneralOptProblem::SetFunctionMin(uint index, vector<double> minPoint, double minValue)
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  mFunctions[index].mMinimumPoint = minPoint;
  mFunctions[index].mMinimumValue = minValue;
  mFunctions[index].mIsMinimumKnown = true;
}

// ------------------------------------------------------------------------------------------------
vector<double> IGeneralOptProblem::GetMinPoint(uint index) const
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  if (mFunctions[index].mIsMinimumKnown == false)
    throw string("Minimum point is unknown");
  return mFunctions[index].mMinimumPoint;
}

// ------------------------------------------------------------------------------------------------
double IGeneralOptProblem::GetMinValue(uint index) const
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  if (mFunctions[index].mIsMinimumKnown == false)
    throw string("Minimum point is unknown");
  return mFunctions[index].mMinimumValue;
}

// ------------------------------------------------------------------------------------------------
void IGeneralOptProblem::SetFunctionMax(uint index, vector<double> maxPoint, double maxValue)
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  mFunctions[index].mMaximumPoint = maxPoint;
  mFunctions[index].mMaximumValue = maxValue;
  mFunctions[index].mIsMaximumKnown = true;
}

// ------------------------------------------------------------------------------------------------
vector<double> IGeneralOptProblem::GetMaxPoint(uint index) const
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  if (mFunctions[index].mIsMaximumKnown == false)
    throw string("Maximum point is unknown");
  return mFunctions[index].mMaximumPoint;
}

// ------------------------------------------------------------------------------------------------
double IGeneralOptProblem::GetMaxValue(uint index) const
{
  if (index >= mFunctions.size())
    throw string("function index is out of range");
  if (mFunctions[index].mIsMaximumKnown == false)
    throw string("Maximum point is unknown");
  return mFunctions[index].mMaximumValue;
}

// ------------------------------------------------------------------------------------------------
vector<double> IGeneralOptProblem::GetOptimumPoint() const
{
  return mOptimumPoint;
}

// ------------------------------------------------------------------------------------------------
double IGeneralOptProblem::GetOptimumValue() const
{
  return mOptimumValue;
}

// ------------------------------------------------------------------------------------------------
vector<double> IGeneralOptProblem::ComputeDerivatives(int index, const vector<double>& y) const
{
  throw string("Derivatives is unknown");
}

// ------------------------------------------------------------------------------------------------
IGeneralOptProblem::IGeneralOptProblem()
{
  mDimension = 0;
  mProblemIndex = -1;
  mFunctionNumber = 0;
}

// ------------------------------------------------------------------------------------------------
IGeneralOptProblem::IGeneralOptProblem(int dim, vector<double> loBound, vector<double> upBound,
  int probIndex) :
  mDimension(dim), mLoBound(loBound), mUpBound(upBound), mProblemIndex(probIndex),
  mFunctionNumber(0)
{ }

// ------------------------------------------------------------------------------------------------
int IGeneralOptProblem::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
void IGeneralOptProblem::GetBounds(vector<double>& lb, vector<double>& ub) const
{
  lb = mLoBound;
  ub = mUpBound;
}