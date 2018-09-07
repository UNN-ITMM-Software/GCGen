#pragma once

#include <vector>
#include <string>

using std::vector;
using std::string;

using uint = unsigned;

enum EOptFunctionParameter { ofpLipschitz, ofpMinimum, ofpMaximum, ofpDerivatives };

/// Базовый класс для задач оптимизации
class IGeneralOptProblem
{
protected:
  // описание отдельной функции в задаче оптимизации
  struct TOptFunction
  {
    // размерность функции
    int mDimension;
    // координаты точки оптимума (глобального минимума)
    vector<double> mMinimumPoint;
    // значение в точке оптимума
    double mMinimumValue;
    // координаты точки глобального максимума
    vector<double> mMaximumPoint;
    // значение в точке глобального максимума
    double mMaximumValue;
    // значение константы Липшица
    double mLipschitzConstant;

    // заданы координаты и значение глобального минимума
    bool mIsMinimumKnown;
    // заданы координаты и значение глобального максимума
    bool mIsMaximumKnown;
    // задано значение константы Липшица
    bool mIsLipschitzConstantKnown;
    // задано вычисление производной функции
    bool mIsDerivativesKnown;
  };
  // количество функций
  int mFunctionNumber;
  // функции задачи оптимизации
  vector<TOptFunction> mFunctions;
  // индексы критериев в списке функций
  vector<int> mCriterionIndeces;
  // индексы ограничений в списке функций
  vector<int> mConstraintIndeces;

  // размерность задачи оптимизации
  int mDimension;
  // Область определения задачи оптимизации в виде гиперинтервала
  // (mLoBound - левый нижний угол, mUpBound - правый верхний угол)
  vector<double> mLoBound;
  vector<double> mUpBound;

  // координаты точки оптимума задачи оптимизации (глобального минимума)
  vector<double> mOptimumPoint;
  // значение в точке оптимума в задаче оптимизации
  double mOptimumValue;

  // номер задачи в семействе, -1 - если не принадлежит семейству
  int mProblemIndex;

  /// Выяснить, что задано для функции с индексом index
  bool GetStatus(uint index, enum EOptFunctionParameter param) const;

  /// Задать значение константы Липшица функции с индексом index
  void SetLipschitzConstant(uint index, double lipConst);
  /// Вернуть значение константы Липшица функции с индексом index
  /// Если константа Липшица не задана, бросается исключение
  double GetLipschitzConstant(uint index) const;

  /// Задать координаты и значение глобального минимума функции с индексом index
  void SetFunctionMin(uint index, vector<double> minPoint, double minValue);
  /// Вернуть координаты глобального минимума функции с индексом index
  /// Если минимальное значение функции не задано, бросается исключение
  vector<double> GetMinPoint(uint index) const;
  /// Вернуть значение глобального минимума функции с индексом index
  /// Если минимальное значение функции не задано, бросается исключение
  double GetMinValue(uint index) const;
  /// Задать координаты и значение глобального максимума функции с индексом index
  void SetFunctionMax(uint index, vector<double> maxPoint, double maxValue);
  /// Вернуть координаты глобального максимума функции с индексом index
  /// Если максимальное значение функции не задано, бросается исключение
  vector<double> GetMaxPoint(uint index) const;
  /// Вернуть значение глобального максимума функции с индексом index
  /// Если максимальное значение функции не задано, бросается исключение
  double GetMaxValue(uint index) const;

  /// Вернуть координаты глобального минимума
  vector<double> GetOptimumPoint() const;
  /// Вернуть значение глобального минимума
  double GetOptimumValue() const;

  /// Вычислить значение функции с индексом index в точке y
  virtual double Compute(int index, const vector<double>& y) const = 0;

  /// Вычислить производные функции с индексом index в точке y
  /// Базовая версия бросает исключение, в потомке м.б. переопределена
  ///   при этом должен быть установлен в true флаг mIsDerivativesKnown
  virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;

  IGeneralOptProblem();
public:
  IGeneralOptProblem(int dim, vector<double> loBound, vector<double> upBound, int probIndex = -1);
  /// Вернуть размерность задачи
  int GetDimension() const;
  /// Вернуть область определения функционала в виде гиперинтервала
  /// (lb - левый нижний угол, ub - правый верхний угол)
  void GetBounds(vector<double>& lb, vector<double>& ub) const;
};