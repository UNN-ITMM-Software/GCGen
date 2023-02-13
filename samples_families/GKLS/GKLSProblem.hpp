#pragma once

/******************************************************************************/
/*        GKLS-Generator of Classes of ND  (non-differentiable),              */
/*                                 D  (continuously differentiable), and      */
/*                                 D2 (twice continuously differentiable)     */
/*                     Test Functions for Global Optimization                 */
/*                                                                            */
/*   Authors:                                                                 */
/*                                                                            */
/*   M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev                         */
/*                                                                            */
/*   (C) 2002-2005                                                            */
/*                                                                            */
/*	 References:                                                              */
/*                                                                            */
/*   1. M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev (2003),              */
/*   Algorithm 829: Software for Generation of Classes of Test Functions      */
/*   with Known Local and Global Minima for Global Optimization.              */
/*   ACM Transactions on Mathematical Software, Vol. 29, no. 4, pp. 469-480.  */
/*                                                                            */
/*   2. D.Knuth (1997), The Art of Computer Programming, Vol. 2:              */
/*   Seminumerical Algorithms (Third Edition). Reading, Massachusetts:        */
/*   Addison-Wesley.                                                          */
/*                                                                            */
/*   The software constructs a convex quadratic function (paraboloid) and then*/
/*   systematically distorts randomly selected parts of this function         */
/*   by polynomials in order to introduce local minima and to construct test  */
/*   functions which are non-differentiable (ND-type), continuously           */
/*   differentiable (D-type), and twice continuously differentiable (D2-type) */
/*   at the feasible region.                                                  */
/*                                                                            */
/*   Each test class is defined by the following parameters:                  */
/*  (1) problem dimension                                                     */
/*  (2) number of local minima including the paraboloid min and the global min*/
/*  (3) global minimum value                                                  */
/*  (3) distance from the paraboloid vertex to the global minimizer           */
/*  (4) radius of the attraction region of the global minimizer               */
/******************************************************************************/


#include "IOptProblem.hpp"
#include "GKLSOption.hpp"

class TGKLSProblem : public IOptProblem
{
protected:
  /// Parameters for curent problem
  GKLSOption option;


  /// Set the parameters
  void SetOptions();


  /// Compute the value of the function number at the point y
  virtual double Compute(int index, const vector<double>& y) const;
  /// Compute the derivatives of the function number index  at the point y
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;

  double CalculateNDFunction(const double* x) const;
  double CalculateDFunction(const double* x) const;
  double CalculateD2Function(const double* x) const;
  double CalculateDFunctionDeriv(unsigned var_j, const double* x) const;
  double CalculateD2FunctionDeriv1(unsigned var_j, const double* x) const;
  double CalculateD2FunctionDeriv2(unsigned var_j, unsigned var_k, const double* x) const;
  int CalculateDFunctionGradient(const double* x, double* g) const;
  int CalculateD2FunctionGradient(const double* x, double* g) const;
  int CalculateD2FunctionHessian(const double* x, double** h) const;

public:
  /// The index is in the range [1,100]
  TGKLSProblem(int problemIndex = 1, int dim = 2, GKLSClass type = Simple, GKLSFuncionType functionType = TD);

  ~TGKLSProblem();
};
