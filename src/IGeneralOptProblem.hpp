#pragma once

#include <vector>
#include <string>

using std::vector;
using std::string;

using uint = unsigned;

enum EOptFunctionParameter { ofpLipschitz, ofpMinimum, ofpMaximum, ofpDerivatives };

/// Base class for optimization problems
class IGeneralOptProblem
{
protected:
  // description of a function in the optimization problem
  struct TOptFunction
  {
    // Dimension
    int mDimension;
    // Minimizer
    vector<double> mMinimumPoint;
    // Minimum value
    double mMinimumValue;
    // Maximizer
    vector<double> mMaximumPoint;
    // Maximum value
    double mMaximumValue;
    // Lipschitz constant
    double mLipschitzConstant;

    // The minimizer and the global minimum value are given
    bool mIsMinimumKnown;
	// The maximizer and the global maximum value are given
    bool mIsMaximumKnown;
    // The value of the Lipschitz constant is given
    bool mIsLipschitzConstantKnown;
    // The function derivative is known
    bool mIsDerivativesKnown;
  };
  // number of functions
  int mFunctionNumber;
  // functions
  vector<TOptFunction> mFunctions;
  // indexes of criteria in the list of functions
  vector<int> mCriterionIndeces;
  // indexes of constraints in the list of functions
  vector<int> mConstraintIndeces;

  // problem dimension
  int mDimension;
  // The problem domain in the form of a hiperinterval 
  // (mLoBound - lower left corner, mUpBound - upper right corner)
  vector<double> mLoBound;
  vector<double> mUpBound;

  // global minimizer
  vector<double> mOptimumPoint;
  // global minimum value
  double mOptimumValue;

  // number of problem in the series, -1 for a single problem
  int mProblemIndex;

  /// Выяснить, что задано для функции с индексом index
  bool GetStatus(uint index, enum EOptFunctionParameter param) const;

  /// What is given for function number index
  void SetLipschitzConstant(uint index, double lipConst);
  /// Lipschitz constant for the function number index
  /// If not specified, throws an exception
  double GetLipschitzConstant(uint index) const;

  /// Set the global minimazer and the global minimum value for function number index
  void SetFunctionMin(uint index, vector<double> minPoint, double minValue);
  /// Get the global minimizer for function number index
  /// If not specified, throws an exception
  vector<double> GetMinPoint(uint index) const;
  /// Get the global minimum value for function number index
  /// If not specified, throws an exception
  double GetMinValue(uint index) const;
  /// Set the global maximazer and the global maximum value for function number index
  void SetFunctionMax(uint index, vector<double> maxPoint, double maxValue);
  /// Get the global maximizer for function number index
  /// If not specified, throws an exception
  vector<double> GetMaxPoint(uint index) const;
  /// Get the global maximum value for function number index
  /// If not specified, throws an exception
  double GetMaxValue(uint index) const;

  /// Get global minimizer
  vector<double> GetOptimumPoint() const;
  /// Get global minimum value
  double GetOptimumValue() const;

  /// Compute the value of the function number index at the point 
  virtual double Compute(int index, const vector<double>& y) const = 0;

  /// Compute the derivatives of the function number index at the point y
  /// If not specified, throws an exception
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;

  IGeneralOptProblem();
public:
  IGeneralOptProblem(int dim, vector<double> loBound, vector<double> upBound, int probIndex = -1);
  /// Get dimension
  int GetDimension() const;
  /// Get domain 
  // (lb - lower left corner, ub - upper right corner)
  void GetBounds(vector<double>& lb, vector<double>& ub) const;
};