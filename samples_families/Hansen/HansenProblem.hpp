#pragma once

#include "IOptProblem.hpp"
#include <cmath>

class THansenProblem0 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 0.0;
    res = pow(y[0], 6) / 6.0 - 52.0 / 25.0*pow(y[0], 5) + 39.0 / 80.0*pow(y[0], 4) +
      71.0 / 10.0*pow(y[0], 3) - 79.0 / 20.0*pow(y[0], 2) - y[0] + 0.1;
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 0.0;
    res = pow(y[0], 5) - 52.0 / 5.0 * pow(y[0], 4) + 39.0 / 20.0 * pow(y[0], 3) +
      213.0 / 10.0 * pow(y[0], 2) - 79.0 / 10.0 * y[0] - 1.0;
    return{ res };
  }
public:
  THansenProblem0() : IOptProblem()
  {
    mDimension = 1;
    mLoBound = { -1.5 };
    mUpBound = { 11.0 };
    mProblemIndex = 0;
    mOptimumPoint = { 10 };
    mOptimumValue = -29763.23;
    TOptFunction func;
    func.mDimension = mDimension;
    func.mMinimumPoint = mOptimumPoint;
    func.mMinimumValue = mOptimumValue;
    func.mIsMinimumKnown = true;
    func.mIsDerivativesKnown = true;
    mFunctions.push_back(func);
    mFunctionNumber = 1;
    mFunctionIndex = 0;
    SetLipschitzConstant(13869.4);
    SetFunctionMax({ 1.41421 }, 2.3847);
  }
};

class THansenProblem1 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 0.0;
    res = sin(y[0]) + sin(10.0 * y[0] / 3.0);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 0.0;
    res = cos(y[0]) + 10.0 / 3.0 * cos(10.0 * y[0] / 3.0);
    return{ res };
  }
public:
  THansenProblem1() : IOptProblem(1, { 2.7 }, { 7.5 }, { 5.14574 }, -1.8996, 1)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(3.95245);
    SetFunctionMax({ 6.21731 }, 0.888315);
  }
};

class THansenProblem2 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 0.0;
    for (int i = 1; i <= 5; i++) res = res + i * sin((i + 1)*y[0] + i);
    return -res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 0.0;
    for (int i = 1; i <= 5; i++) res += i * (i + 1) * cos((i + 1)*y[0] + i);
    return{ -res };
  }
public:
  THansenProblem2() : IOptProblem(1, { -10.0 }, { 10.0 }, { -0.49139 }, -12.0312, 2)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(68.4194);
    SetFunctionMax({ -1.1141 }, 14.838);
  }
};

class THansenProblem3 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = (-16.0 * y[0] * y[0] + 24.0 * y[0] - 5.0) * exp(-y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = (16.0 * y[0] * y[0] - 56.0 * y[0] + 29.0) * exp(-y[0]);
    return{ res };
  }
public:
  THansenProblem3() : IOptProblem(1, { 1.9 }, { 3.9 }, { 2.86803 }, -3.85045, 3)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(2.93753);
    SetFunctionMax({ 1.9 }, -2.5666);
  }
};

class THansenProblem4 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = -(-3.0 * y[0] + 1.4) * sin(18.0 * y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 3.0 * sin(18.0 * y[0]) + 18.0 * (3.0 * y[0] - 1.4) * cos(18.0 * y[0]);
    return{ res };
  }
public:
  THansenProblem4() : IOptProblem(1, { 0.0 }, { 1.2 }, { 0.96609 }, -1.48907, 4)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(35.4653);
    SetFunctionMax({ 1.13904 }, 2.01028);
  }
};

class THansenProblem5 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = -(y[0] + sin(y[0])) * exp(-y[0] * y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = (-1.0 - cos(y[0]) + 2.0 * y[0] * y[0] + 2.0 * y[0] * sin(y[0])) * exp(-y[0] * y[0]);
    return{ res };
  }
public:
  THansenProblem5() : IOptProblem(1, { -10.0 }, { 10.0 }, { 0.67958 }, -0.824239, 5)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(2.0);
    SetFunctionMax({ -0.67958 }, 0.824239);
  }
};

class THansenProblem6 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = sin(y[0]) + sin(10.0 * y[0] / 3.0) + log(y[0]) - 0.84*y[0] + 3.0;
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = cos(y[0]) + 10.0 / 3.0 * cos(10.0 * y[0] / 3.0) + 1.0 / y[0] - 0.84;
    return{ res };
  }
public:
  THansenProblem6() : IOptProblem(1, { 2.7 }, { 7.5 }, { 5.19978 }, -1.60131, 6)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(4.44012);
    SetFunctionMax({ 2.7 }, 2.56475);
  }
};

class THansenProblem7 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 0;
    for (int i = 1; i <= 5; i++)
      res += i*cos((i + 1)*y[0] + i);
    return -res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 0.0;
    for (int i = 1; i <= 5; i++)
      res += i * (i + 1) * sin((i + 1)*y[0] + i);
    return{ res };
  }
public:
  THansenProblem7() : IOptProblem(1, { -10.0 }, { 10.0 }, { -0.80032 }, -14.508, 7)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(69.4801);
    SetFunctionMax({ -1.42513 }, 12.8709);
  }
};

class THansenProblem8 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = sin(y[0]) + sin(2.0 / 3.0 * y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = cos(y[0]) + 2.0 / 3.0 * cos(2.0 / 3.0 * y[0]);
    return{ res };
  }
public:
  THansenProblem8() : IOptProblem(1, { 3.1 }, { 20.4 }, { 17.0392 }, -1.90596, 8)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(1.66667);
    SetFunctionMax({ 20.4 }, 1.85895);
  }
};

class THansenProblem9 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = -y[0] * sin(y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = -sin(y[0]) - y[0] * cos(y[0]);
    return{ res };
  }
public:
  THansenProblem9() : IOptProblem(1, { 0.0 }, { 10.0 }, { 7.97867 }, -7.91673, 9)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(9.63171);
    SetFunctionMax({ 10.0 }, 5.44021);
  }
};

class THansenProblem10 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 2.0 * cos(y[0]) + cos(2 * y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = -2.0 * sin(y[0]) - 2.0 * sin(2.0 * y[0]);
    return{ res };
  }
public:
  THansenProblem10() : IOptProblem(1, { -1.57 }, { 6.28 }, { 4.18879 }, -1.5, 10)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(3.52035);
    SetFunctionMax({ 0.0 }, 3.0);
  }
};

class THansenProblem11 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = pow(sin(y[0]), 3) + pow(cos(y[0]), 3);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 3.0 * sin(y[0]) * cos(y[0]) * (sin(y[0]) - cos(y[0]));
    return{ res };
  }
public:
  THansenProblem11() : IOptProblem(1, { 0.0 }, { 6.28 }, { 4.71239 }, -1.0, 11)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(2.12132);
    SetFunctionMax({ 0.0 }, 1.0);
  }
};

class THansenProblem12 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double sgn, res;
    sgn = (y[0] * y[0] - 1.0 < 0.0) ? -1.0 : 1.0;
    res = -pow(y[0] * y[0], 1.0 / 3.0) + sgn * pow(sgn * (y[0] * y[0] - 1.0), 1.0 / 3.0);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double sgn, res;
    sgn = (y[0] * y[0] - 1.0 < 0.0) ? -1.0 : 1.0;
    res = -2.0 / 3.0 * pow(y[0], -1.0 / 3.0) + 1.0 / 3.0 * pow(sgn * (y[0] * y[0] - 1.0), -2.0 / 3.0) * 2.0 * y[0];
    return{ res };
  }
public:
  THansenProblem12() : IOptProblem(1, { 0.001 }, { 0.99 }, { 0.70711 }, -1.5874, 12)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(8.5);
    SetFunctionMax({ 0.001 }, -1.01);
  }
};

class THansenProblem13 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = -exp(-y[0]) * sin(2.0 * acos(-1.0) * y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = -exp(-y[0]) * (2.0 * acos(-1.0) * cos(2.0 * acos(-1.0) * y[0]) - sin(2.0 * acos(-1.0) * y[0]));
    return{ res };
  }
public:
  THansenProblem13() : IOptProblem(1, { 0.0 }, { 4.0 }, { 0.22488 }, -0.788685, 13)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(6.28319);
    SetFunctionMax({ 0.72488 }, 0.478362);
  }
};

class THansenProblem14 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = (y[0] * y[0] - 5.0 * y[0] + 6.0) / (y[0] * y[0] + 1);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = (5.0 * y[0] * y[0] - 10.0 * y[0] - 5.0) / (y[0] * y[0] + 1.0) / (y[0] * y[0] + 1.0);
    return{ res };
  }
public:
  THansenProblem14() : IOptProblem(1, { -5.0 }, { 5.0 }, { 2.41421 }, -0.0355339, 14)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(6.3726);
    SetFunctionMax({ -0.41421 }, 7.03553);
  }
};

class THansenProblem15 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = 2.0 * (y[0] - 3.0)*(y[0] - 3.0) + exp(y[0] * y[0] / 2.0);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 4.0 * (y[0] - 3.0) + y[0] * exp(y[0] * y[0] / 2.0);
    return{ res };
  }
public:
  THansenProblem15() : IOptProblem(1, { -3.0 }, { 3.0 }, { 1.59072 }, 7.51592, 15)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(294.051);
    SetFunctionMax({ -3.0 }, 162.017);
  }
};

class THansenProblem16 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = pow(y[0], 6) - 15.0 * pow(y[0], 4) + 27.0 * y[0] * y[0] + 250.0;
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = 6.0 * pow(y[0], 5) - 60.0 * pow(y[0], 3) + 54.0 * y[0];
    return{ res };
  }
public:
  THansenProblem16() : IOptProblem(1, { -4.0 }, { 4.0 }, { 3.0 }, 7.0, 16)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(2520.0);
    SetFunctionMax({ 4.0 }, 938.0);
  }
};

class THansenProblem17 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = (y[0] <= 3.0) ? (y[0] - 2.0)*(y[0] - 2.0) : 2.0 * log(y[0] - 2.0) + 1.0;
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = (y[0] <= 3.0) ? 2.0 * (y[0] - 2.0) : 2.0 / (y[0] - 2.0);
    return{ res };
  }
public:
  THansenProblem17() : IOptProblem(1, { 0.0 }, { 6.0 }, { 2.0 }, 0, 17)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(4.0);
    SetFunctionMax({ 0.0 }, 4.0);
  }
};

class THansenProblem18 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = -y[0] + sin(3.0 * y[0]) - 1.0;
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = -1.0 + 3.0 * cos(3.0 * y[0]);
    return{ res };
  }
public:
  THansenProblem18() : IOptProblem(1, { 0.0 }, { 6.5 }, { 5.87287 }, -7.81567, 18)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(4.0);
    SetFunctionMax({ 0.41032 }, -0.467511);
  }
};

class THansenProblem19 : public IOptProblem
{
  virtual double Compute(int index, const vector<double>& y) const
  {
    double res = -(y[0] - sin(y[0])) * exp(-y[0] * y[0]);
    return res;
  }
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const
  {
    double res = (-1.0 + cos(y[0]) + 2.0 * y[0] * y[0] - 2.0 * y[0] * sin(y[0])) * exp(-y[0] * y[0]);
    return{ res };
  }
public:
  THansenProblem19() : IOptProblem(1, { -10.0 }, { 10.0 }, { 1.19514 }, -0.0634905, 19)
  {
    mFunctions[mFunctionIndex].mIsDerivativesKnown = true;
    SetLipschitzConstant(0.0962709);
    SetFunctionMax({ -1.19514 }, 0.0634905);
  }
};
