#pragma once

#include "IGeneralOptProblem.hpp"

class IOptProblem : public IGeneralOptProblem
{
protected:
  /// Сriterion number in the vector mFunctions
  uint mFunctionIndex;
  /// Set the value of Lipschitz constant
  void SetLipschitzConstant(double lipConst);
  /// Set global maximizer and global maximum value
  void SetFunctionMax(vector<double> maxPoint, double maxValue);
  IOptProblem();
public:
  IOptProblem(int dim, vector<double> loBound, vector<double> upBound,
    vector<double> optPoint, double optVal, int probIndex = -1);

  /// Get global minimizer
  vector<double> GetOptimumPoint() const;
  /// Get global minimum value
  double GetOptimumValue() const;

  /// What is specified for the objective function
  bool GetStatus(enum EOptFunctionParameter param) const;
  /// Get global maximizer
  vector<double> GetMaxPoint() const;
  /// Get global maximum value
  double GetMaxValue() const;
  /// Get the value of Lipschitz constant
  double GetLipschitzConstant() const;

  /// Compute the value of the objective function at the point y
  double ComputeFunction(const vector<double>& y) const;

  /// Compute the value of the objective function derivatives at the point y
  vector<double> ComputeFunctionDerivatives(const vector<double>& y) const;
};
