#pragma once

/* Two-dimensional multiextremal test function of Vladimir A. Grishagin

   See:  Grishagin, V.A., Operating Characteristics of Some Global Search Algorithms.
         in: Problems in Random Search, Riga: Zinatne, 1978, issue 7, pp. 198--206;
   also: Strongin, R.G., and Sergeyev, Ya.D., Global Optimization with Non-Convex
         Constraints: Sequential and Parallel Algorithms. Dordrecht: Kluwer, 2000.
*/

#include "IOptProblem.hpp"
#include "GrishaginOption.hpp"

class TGrishaginProblem : public IOptProblem
{
protected:
  /// Current problem parameters
  GrishaginOption option;

  double rndm20(unsigned char k[]);
  void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);
  /// Set parameters
  void SetOptions();

  /// Compute x-derivaive
  double CalculateXDerivative(const vector<double>& y) const;
  /// Compute y-derivaive
  double CalculateYDerivative(const vector<double>& y) const;

  /// Compute the function number index at the point y
  virtual double Compute(int index, const vector<double>& y) const;
  /// Compute the derivatives of the function number index at the point y
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;

public:
  /// ProblemIndex is in the range [1,100]
  TGrishaginProblem(int problemIndex = 1);
};
