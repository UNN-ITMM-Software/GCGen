#pragma once

#include "IConstrainedOptProblem.hpp"
#include "GrishaginOption.hpp"

// Пример класса с определением задачи условной оптимизации
class GrishaginConstrainedProblem : public IConstrainedOptProblem
{
protected:

  /// Параметры для текущей задачи
  vector<GrishaginOption> options;

  double rndm20(unsigned char k[]);
  void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);
  /// Задать параметры функции
  void SetOptions(int index);

  /// Вычисление производной по нулевой координате
  double CalculateXDerivative(int index, const vector<double>& y) const;
  /// Вычисление производной по первой координате
  double CalculateYDerivative(int index, const vector<double>& y) const;

  /// Вычислить значение функции с индексом index в точке y
  virtual double Compute(int index, const vector<double>& y) const;

  /// Вычислить производные функции с индексом index в точке y
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;


public:
  GrishaginConstrainedProblem(EConstrainedProblemType problemType = cptInFeasibleDomain,
    double fraction = 0.5, int activeConstrNum = 0, int problemIndex = 1);
};