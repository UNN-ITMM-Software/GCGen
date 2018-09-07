#pragma once

#include "IGeneralOptProblem.hpp"


// Способы формирования задачи с ограничениями
enum EConstrainedProblemType {
  cptNormal, cptInFeasibleDomain, cptOutFeasibleDomain,
  cptOnFeasibleBorder
};

// Схемы вычисления ограничений
enum EConstraintComputationType { cctAllConstraints, cctIndexScheme };

class IConstrainedOptProblem : public IGeneralOptProblem
{
protected:
  // тип формирования задачи с ограничениями
  EConstrainedProblemType mProblemType;
  // доля допустимой области
  double mFeasibleDomainFraction;
  // количество активных ограничений (только для варианта OnFeasibleBorder)
  // по умолчанию = 0, т.е. параметр не установлен
  int mActiveConstraintNumber;

  /// Сдвиг ограничений (RHS)
  std::vector<double> mQ;
  /// Коэфициенты масштабирования для ограничений
  std::vector<double> mZoomRatios;
  /// Покоординатный сдвиг ограничений
  std::vector<std::vector<double> > mShift;
  /// Покоординатный сдвиг глобального минимума на границу
  std::vector<double> mBoundaryShift;
  /// Коэффициенты изменения
  std::vector<double> mImprovementCoefficients;
  /// Точность поиска ближайшей точки на границе области (степень 0.5)
  int mBoundarySearchPrecision;


  /// Нужно ли масштабировать задачу
  bool mIsZoom;
  /// Нужно ли сдвигать функции ограничений
  bool mIsShift;
  /// Нужно ли сдвигать оптимум целевой функции на границу
  bool mIsBoundaryShift;
  /// Одно общее дельта(доля области) для всех функций или для каждой функции свое
  bool mIsTotalDelta;
  /// Изменять ли целевую функцию путем прибаления функционала от ограничений
  bool mIsImprovementOfTheObjective;

  void InitProblem();

  // выполнить преобразование значения целевой функции,
  //   исходя из значения mProblemType и параметров задачи
  double TransformValue(double val, vector<double> point, int index = 0) const;
  // выполнить преобразование координат точки целевой функции,
  //   исходя из значения mProblemType и параметров задачи
  vector<double> TransformPoint(vector<double> point, int index = 0) const;
  // выполнить преобразование значения ограничения,
  //   исходя из значения mProblemType и параметров задачи
  double TransformConstraintValue(double val, int index) const;
  // выполнить преобразование координат точки ограничения,
  //   исходя из значения mProblemType и параметров задачи
  vector<double> TransformConstraintPoint(vector<double> point, int index) const;

  /// Максимум значений всех ограничений в точке
  double MaxFunctionCalculate(std::vector<double> y);

  /** Вычисляет сдвиг ограничений по доле допустимой области
  \param[in] delta - заданная доля допустимой области
  \param[in] m - количество шагов сетки при определении доли допустимой области
  \param[in] Epsilon - заданная точность вычисления доли допустимой области
  \param[in] maxM - максимальное число испытаний
  \return сдвиг ограничений
  */
  double CalculateRHS(double delta, int m = 100, double Epsilon = 0.01, int maxM = 10000000);

  /// Задает коэфициенты масштабирования для функций ограничений
  virtual void SetZoom();
  /// Задает сдвиг к глобальному минимуму для функций ограничений
  virtual void SetShift();
  /// Задает сдвиг глобального минимума целевой функции на границу области
  virtual vector<double>  SetBoundaryShift();

public:
  IConstrainedOptProblem(int dim, vector<double> loBound, vector<double> upBound,
    int probIndex = -1, EConstrainedProblemType problemType = cptNormal,
    double fraction = 1, int activeConstrNum = 0);

  /// Вернуть координаты глобального минимума (в допустимой области)
  vector<double> GetOptimumPoint() const;
  /// Вернуть значение глобального минимума (в допустимой области)
  double GetOptimumValue() const;

  /// Выяснить, что задано для целевой функции
  bool GetStatus(enum EOptFunctionParameter param) const;
  /// Вернуть координаты глобального максимума целевой функции
  vector<double> GetMaxPoint() const;
  /// Вернуть значение глобального максимума целевой функции
  double GetMaxValue() const;
  /// Вернуть значение константы Липшица целевой функции
  double GetLipschitzConstant() const;

  /// Вычислить значение целевой функции в точке y
  double ComputeFunction(const vector<double>& y) const;
  /// Вычислить производные целевой функции в точке y
  /// Базовая версия бросает исключение, в потомке м.б. переопределена
  ///   при этом должен быть установлен в true флаг mIsDerivativesKnown
  vector<double> ComputeFunctionDerivatives(const vector<double>& y) const;

  /// Вернуть число ограничений
  int GetConstraintsNumber() const;

  /// Выяснить, что задано для ограничения
  bool GetConstraintStatus(int index, enum EOptFunctionParameter param) const;
  /// Вернуть значение константы Липшица ограничения
  double GetConstraintLipschitzConstant(int index) const;

  /// Вычислить значение ограничения в точке y
  double ComputeConstraint(int index, const vector<double>& y) const;

  /// Вычислить значения ограничений в точке y, index - номер последнего вычисленного ограничения
  vector<double> ComputeConstraints(const vector<double>& y,
    EConstraintComputationType t, int &index) const;

  /// Вычислить производные ограничения в точке y
  vector<double> ComputeConstraintDerivatives(int index, const vector<double>& y) const;

  virtual ~IConstrainedOptProblem();
};