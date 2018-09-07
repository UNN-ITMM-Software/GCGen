#pragma once

#include "IGeneralOptProblem.hpp"

class IOptProblem : public IGeneralOptProblem
{
protected:
  /// Номер критерия в векторе mFunctions
  uint mFunctionIndex;
  /// Задать значение константы Липшица
  void SetLipschitzConstant(double lipConst);
  /// Задать координаты и значение глобального максимума
  void SetFunctionMax(vector<double> maxPoint, double maxValue);
  IOptProblem();
public:
  IOptProblem(int dim, vector<double> loBound, vector<double> upBound,
    vector<double> optPoint, double optVal, int probIndex = -1);

  /// Вернуть координаты глобального минимума
  vector<double> GetOptimumPoint() const;
  /// Вернуть значение глобального минимума
  double GetOptimumValue() const;

  /// Выяснить, что задано для целевой функции
  bool GetStatus(enum EOptFunctionParameter param) const;
  /// Вернуть координаты глобального максимума
  vector<double> GetMaxPoint() const;
  /// Вернуть значение глобального максимума
  double GetMaxValue(int index) const;
  /// Вернуть значение константы Липшица
  double GetLipschitzConstant() const;

  /// Вычислить значение целевой функции в точке y
  double ComputeFunction(const vector<double>& y) const;

  /// Вычислить производные целевой функции в точке y
  vector<double> ComputeFunctionDerivatives(const vector<double>& y) const;
};
