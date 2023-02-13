#pragma once

#include "IConstrainedOptProblem.hpp"
#include "GKLSOption.hpp"

// Example of class for constrained problem
class TGKLSConstrainedProblem : public IConstrainedOptProblem
{
protected:

  /// Parameters
  vector<GKLSOption> options;

  /// Set parameters
  void SetOptions(int problemIndex);

  /// Compute the function number index at the point y
  virtual double Compute(int problemIndex, const vector<double>& y) const;

  /// Compute the derivatives of the function number index at the point y
  virtual vector<double> ComputeDerivatives(int problemIndex, const vector<double>& y) const;


  double CalculateNDFunction(int problemIndex, const double* x) const;
  double CalculateDFunction(int problemIndex, const double* x) const;
  double CalculateD2Function(int problemIndex, const double* x) const;
  double CalculateDFunctionDeriv(int problemIndex, unsigned var_j, const double* x) const;
  double CalculateD2FunctionDeriv1(int problemIndex, unsigned var_j, const double* x) const;
  double CalculateD2FunctionDeriv2(int problemIndex, unsigned var_j, unsigned var_k, const double* x) const;
  int CalculateDFunctionGradient(int problemIndex, const double* x, double* g) const;
  int CalculateD2FunctionGradient(int problemIndex, const double* x, double* g) const;
  int CalculateD2FunctionHessian(int problemIndex, const double* x, double** h) const;


public:
  TGKLSConstrainedProblem(EConstrainedProblemType problemType = cptInFeasibleDomain,
    double fraction = 0.5, int activeConstrNum = 0, int problemIndex = 1, int dim = 2,
    GKLSClass type = Simple, GKLSFuncionType functionType = TD);

  ~TGKLSConstrainedProblem();
};