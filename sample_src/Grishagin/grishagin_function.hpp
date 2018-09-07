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
  /// Параметры для текущей задачи
  GrishaginOption option;

  double rndm20(unsigned char k[]);
  void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);
  /// Задать параметры функции
  void SetOptions();

  /// Вычисление производной по нулевой координате
  double CalculateXDerivative(const vector<double>& y) const;
  /// Вычисление производной по первой координате
  double CalculateYDerivative(const vector<double>& y) const;

  /// Вычислить значение функции с индексом index в точке y
  virtual double Compute(int index, const vector<double>& y) const;
  /// Вычислить производные функции с индексом index в точке y
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;

public:
  /// Принимает индекс функции от 1 до 100 включительно
  TGrishaginProblem(int problemIndex = 1);
};
