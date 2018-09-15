#pragma once

#include "IConstrainedOptProblem.hpp"
#include "GrishaginOption.hpp"

// Example of class for constrained problem
class GrishaginConstrainedProblem : public IConstrainedOptProblem
{
protected:

  /// Parameters
  vector<GrishaginOption> options;

  double rndm20(unsigned char k[]);
  void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);
  /// Set parameters
  void SetOptions(int index);

  /// Compute x-derivaive
  double CalculateXDerivative(int index, const vector<double>& y) const;
  /// Compute y-derivaive
  double CalculateYDerivative(int index, const vector<double>& y) const;

  /// Compute the function number index at the point y
  virtual double Compute(int index, const vector<double>& y) const;

  /// Compute the derivatives of the function number index at the point y
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;


public:
  GrishaginConstrainedProblem(EConstrainedProblemType problemType = cptInFeasibleDomain,
    double fraction = 0.5, int activeConstrNum = 0, int problemIndex = 1);
};