#include "IConstrainedOptProblem.hpp"
#include <algorithm>

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::MaxFunctionCalculate(vector<double> y)
{
  if (mConstraintIndeces.size() > 0)
  {
    vector<double> x = TransformConstraintPoint(y, 0);

    double res = Compute(mConstraintIndeces[0], x);
    for (uint i = 1; i < mConstraintIndeces.size(); i++)
    {
      x = TransformConstraintPoint(y, i);

      double f = Compute(mConstraintIndeces[i], x);
      if (res < f)
        res = f;
    }
    return res;
  }
  return 0;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::CalculateRHS(double delta, int m, double Epsilon,
  int maxM)
{
  double rhs = 0;
  double hmin = mOptimumValue;
  double hmax = hmin;
  double d = 0;

  // многомерная решетка, в узлах - испытания
  int* size = new int[mDimension]; // кол-во точек на размерность
  double* step = new double[mDimension]; // шаг по каждой размерности
  int sumn = 1; // число испытаний

  double* a = mLoBound.data();
  double* b = mUpBound.data();

  double multiplyingLength = 1;
  for (int i = 0; i < mDimension; i++)
  {
    d = (b[i] - a[i]);
    size[i] = (int)ceil(d / Epsilon) + 1;
    step[i] = d / (size[i] - 1);
    multiplyingLength = multiplyingLength * d;
    sumn *= (size[i]);
  }

  if ((sumn > maxM) || (sumn <= 0))
  {
    multiplyingLength = multiplyingLength / maxM;
    Epsilon = pow(multiplyingLength, 1.0 / (double)mDimension);
    sumn = 1;
    multiplyingLength = 1;

    for (int i = 0; i < mDimension; i++)
    {
      d = (b[i] - a[i]);
      size[i] = (int)ceil(d / Epsilon) + 1;
      step[i] = d / (size[i] - 1);
      sumn *= (size[i]);
    }
  }

  double* f = new double[sumn]; // значение функции
  vector<double> yArray(mDimension);

  for (int i = 0; i < sumn; i++)
  {
    double w;
    int z = i;
    // вычисляем координаты точки испытания
    for (int j = 0; j < mDimension; j++)
    {
      w = z % size[j]; // определяем номер узла на i-ой размерности
      yArray[j] = a[j] + w * step[j];//левая граница + номер узла на i-ой размерности * шаг на
                                // i-ой размерности
      z = z / size[j]; // для вычисления номера узла на следующей размерности
    }
    // проводим испытание
    f[i] = MaxFunctionCalculate(yArray);
    hmax = std::max(f[i], hmax);
    hmin = std::min(f[i], hmin);
  }

  double* h1 = new double[m];
  double* h2 = new double[m];
  int* p = new int[m];
  int* s = new int[m];

  double deltah = (hmax - hmin) / m;

  for (int i = 0; i < m; i++)
  {
    h1[i] = hmin + i * deltah;
    h2[i] = hmin + (i + 1) * deltah;
    p[i] = 0;
    s[i] = 0;
  }

  for (int i = 0; i < sumn; i++)
    for (int j = 0; j < m; j++)
      if ((f[i] >= h1[j]) && (f[i] <= h2[j]))
      {
        p[j] ++;
        break;
      }

  s[0] = p[0];
  for (int i = 1; i < m; i++)
  {
    s[i] = s[i - 1] + p[i];
  }

  double smax = s[m - 1];
  double g = delta * smax;
  for (int i = 0; i < m; i++)
  {
    if (s[i] >= g)
    {
      rhs = h2[i];
      break;
    }
  }

  double dm = delta;
  if (dm == 0)
    dm += 0.1;
  dm = dm * (hmax - hmin);

  double criticalValue = MaxFunctionCalculate(mOptimumPoint);

  if (rhs < criticalValue)
  {
//    std::cout << "Q was changed from " << rhs << " to " << criticalValue + dm << "\n";
    rhs = criticalValue + dm;
  }

  delete[] size;
  delete[] step;
  delete[] f;

  delete[] h1;
  delete[] h2;
  delete[] p;
  delete[] s;

  return rhs;
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::SetZoom()
{
  double AccuracyDouble = 0.00000001;
  double* lower = mLoBound.data();
  double* upper = mUpBound.data();
  double* constraintMin = 0;

  double maxDistanceToBoundary = 0;

  for (int k = 0; k < mDimension; k++)
  {
    if (maxDistanceToBoundary < (mOptimumPoint[k] - lower[k]))
      maxDistanceToBoundary = (mOptimumPoint[k] - lower[k]);
    if (maxDistanceToBoundary < (upper[k] - mOptimumPoint[k]))
      maxDistanceToBoundary = (upper[k] - mOptimumPoint[k]);
  }

  if (fabs(maxDistanceToBoundary) < AccuracyDouble)
  {
    mIsZoom = false;
    for (uint j = 0; j < mConstraintIndeces.size(); j++)
    {
      mZoomRatios[j] = 1;
    }

    return;
  }

  for (uint i = 0; i < mConstraintIndeces.size(); i++)
  {
    double minDistanceToBoundary = upper[0] - lower[0];

    constraintMin = mFunctions[mConstraintIndeces[i]].mMinimumPoint.data();

    for (int k = 0; k < mDimension; k++)
    {
      if (minDistanceToBoundary > (constraintMin[k] - lower[k]))
        minDistanceToBoundary = (constraintMin[k] - lower[k]);
      if (minDistanceToBoundary > (upper[k] - constraintMin[k]))
        minDistanceToBoundary = (upper[k] - constraintMin[k]);
    }

    if (fabs(minDistanceToBoundary) < AccuracyDouble)
    {
      mIsZoom = false;
      for (uint j = 0; j < mConstraintIndeces.size(); j++)
      {
        mZoomRatios[j] = 1;
      }

      return;
    }
    else
    {
      mZoomRatios[i] = maxDistanceToBoundary / minDistanceToBoundary;
    }
  }
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::SetShift()
{
  double* constraintMin = 0;

  for (uint i = 0; i < mConstraintIndeces.size(); i++)
  {
    constraintMin = mFunctions[mConstraintIndeces[i]].mMinimumPoint.data();

    for (int k = 0; k < mDimension; k++)
    {
      mShift[i][k] = mOptimumPoint[k] - constraintMin[k] * mZoomRatios[i];
    }
  }
}

// ------------------------------------------------------------------------------------------------
vector<double>  IConstrainedOptProblem::SetBoundaryShift()
{
  double* objectiveMin = mOptimumPoint.data();
  vector<vector<double>> tempPoint(2 * mDimension); // точки на границе по 2*mDimension направлениям
  double* lower = mLoBound.data(); // верхняя граница допустимой области
  double* upper = mUpBound.data(); // нижняя граница допустимой области
  double* outPoint = new double[mDimension], *inPoint = new double[mDimension]; // вспомогательные точки для уточнения точки на границе
  double delta = pow(0.5, 7); // первичный шаг поиска
  vector<double> activeShift(mConstraintIndeces.size());

  for (unsigned j = 0; j < mConstraintIndeces.size(); j++)
  {
    activeShift[j] = 0; // сдвиг ограничений для достижения заданного количества активных ограничений
  }
  for (int j = 0; j < 2 * mDimension; j++)
  {
    tempPoint[j].resize(mDimension);
    for (int k = 0; k < mDimension; k++)
    {
      tempPoint[j][k] = 0;
    }
  }
  for (int j = 0; j < mDimension; j++)
  {
    outPoint[j] = 0;
    inPoint[j] = 0;
  }


  bool isBoundReached = false; // найдена ли точка на границе
  int closestDir = 0, closestDirFull = 0; // координата и направление, в котором была найдена ближайшая точка на границе
  int i = 1;
  unsigned dirNum = 0; // маска: по каким направлениям граничные точки уже найдены

  while ((!isBoundReached) || (dirNum < (unsigned)pow(2, 2 * mDimension) - 1))
  {
    for (int k = 0; k < 2 * mDimension; k++)
    {
      if ((dirNum & (unsigned)pow(2, k)) == (unsigned)pow(2, k)) continue; // пропустить итерацию, если в этом направлении граничная точка найдена
      isBoundReached = false;
      // расчет новой точки
      for (int i = 0; i < mDimension; i++)
      {
        tempPoint[k][i] = objectiveMin[i] + mBoundaryShift[i];// GetOptimumPoint(tempPoint[k]);
      }

      tempPoint[k][k%mDimension] = objectiveMin[k%mDimension] + ((k >= mDimension) ? (-1) : (1)) * delta * i;

      // проверка, не вышли ли за границы области поиска (без учета ограничений)
      if ((tempPoint[k][k % mDimension] - upper[k % mDimension] >= 0) ||
        (tempPoint[k][k % mDimension] - lower[k % mDimension] <= 0))
      {
        if (tempPoint[k][k % mDimension] - upper[k % mDimension] >= 0)
        {
          tempPoint[k][k % mDimension] = upper[k % mDimension];
        }
        else
        {
          tempPoint[k][k % mDimension] = lower[k % mDimension];
        }
        isBoundReached = true;
        closestDir = k % mDimension;
        closestDirFull = k;
        dirNum = dirNum | (unsigned)pow(2, k);
        continue;
      }

      for (uint j = 0; j < mConstraintIndeces.size(); j++)
      {
        if (ComputeConstraint(j, tempPoint[k]) >= 0)
        //if (CalculateFunctionConstrained(tempPoint[k], j) >= 0) // впервые нарушено ограничение в этом направлении
        {
          isBoundReached = true;
          closestDir = k%mDimension;
          closestDirFull = k;

          if (ComputeConstraint(j, tempPoint[k]) == 0)
          //if (CalculateFunctionConstrained(tempPoint[k], j) == 0) // нашли точку строго на границе допустимой области
          {
            break;
          }

          // иначе - вышли за границу -> находим точку ВНУТРИ допустимой области делением пополам
          for (int l = 0; l < mDimension; l++)
          {
            outPoint[l] = tempPoint[k][l];
            inPoint[l] = tempPoint[k][l];
          }
          inPoint[closestDir] = objectiveMin[closestDir] + ((k >= mDimension) ? (-1) : (1)) * delta * (i - 1);

          // точность поиска  - 2^(-mBoundarySearchPrecision)
          while ((abs(inPoint[closestDir] - outPoint[closestDir]) > pow(0.5, mBoundarySearchPrecision)) &&
            (ComputeConstraint(j, tempPoint[k]) != 0))
            //(CalculateFunctionConstrained(tempPoint[k], j) != 0))
          {
            delta = delta / 2;
            tempPoint[k][closestDir] = inPoint[closestDir] + ((k >= mDimension) ? (-1) : (1)) * delta;
            if (ComputeConstraint(j, tempPoint[k]) > 0)
            //if (CalculateFunctionConstrained(tempPoint[k], j) > 0)
            {
              outPoint[closestDir] = tempPoint[k][closestDir];
              tempPoint[k][closestDir] = inPoint[closestDir];
            }
            else if (ComputeConstraint(j, tempPoint[k]) < 0)
            //else if (CalculateFunctionConstrained(tempPoint[k], j) < 0)
            {
              inPoint[closestDir] = tempPoint[k][closestDir];
            }
          }
          break;
        }
      }
      if (isBoundReached)
      {
        dirNum = dirNum | (unsigned)pow(2, k);
        if (mActiveConstraintNumber == 1) // достаточно ближайшей точки с 1 активным ограничением -> не продолжаем перебор направлений
          break;
      }
    }
    if ((isBoundReached) && (mActiveConstraintNumber == 1)) break;// достаточно ближайшей точки с 1 активным ограничением -> не продолжаем увеличение радиуса поиска
    i++;
  }

  double ** deltaForConstraints = new double*[2 * mDimension]; // в каждой найденной точке массив сдвигов каждого из ограничений, чтобы сделать его активным
  double ** sortDeltaForConstraints = new double*[2 * mDimension];
  for (int k = 0; k < 2 * mDimension; k++)
  {
    deltaForConstraints[k] = new double[mConstraintIndeces.size()];
    sortDeltaForConstraints[k] = new double[mConstraintIndeces.size()];
  }
  double diff = 0, minDiff = 50000; // суммарный сдвиг ограничений для достижения заданного количества активных
  if (mActiveConstraintNumber > 1)
  {
    for (int k = 0; k < 2 * mDimension; k++)
    {
      // в каждой найденной точке вычисляем сдвиги ограничений, чтобы сделать их активными
      for (int j = 0; j < mConstraintIndeces.size(); j++)
      {
        if (ComputeConstraint(j, tempPoint[k]) >= -0.03)
        //if (CalculateFunctionConstrained(tempPoint[k], j) >= -0.03)
        {
          deltaForConstraints[k][j] = 0;
          sortDeltaForConstraints[k][j] = 0;
        }
        else
        {
          deltaForConstraints[k][j] = ComputeConstraint(j, tempPoint[k]);
          sortDeltaForConstraints[k][j] = abs(ComputeConstraint(j, tempPoint[k]));

          //deltaForConstraints[k][j] = CalculateFunctionConstrained(tempPoint[k], j);
          //sortDeltaForConstraints[k][j] = abs(CalculateFunctionConstrained(tempPoint[k], j));
        }
      }
      // сортируем сдвиги по возрастанию
      std::sort(sortDeltaForConstraints[k], sortDeltaForConstraints[k] + mConstraintIndeces.size());
      diff = 0;
      // находим суммарный сдвиг ограничений для достижения заданного количества активных
      for (int n = 0; n < mActiveConstraintNumber; n++)
      {
        diff = diff + abs(sortDeltaForConstraints[k][n]);
      }
      // определяем точку с наименьшей величиной суммарного сдвига
      if (diff < minDiff)
      {
        minDiff = diff;
        closestDirFull = k;
      }
    }

    // выставляем сдвиги
    for (int n = 0; n < mActiveConstraintNumber; n++)
    {
      int l = 0;
      while (sortDeltaForConstraints[closestDirFull][n] != abs(deltaForConstraints[closestDirFull][l]))
      {
        l++;
      }
      activeShift[l] = -deltaForConstraints[closestDirFull][l];
    }
  }

  // вычисляем величину сдвига точки глобального минимума
  for (int j = 0; j < mDimension; j++)
  {
    mBoundaryShift[j] = tempPoint[closestDirFull][j] - objectiveMin[j];
  }

  delete[] outPoint;
  delete[] inPoint;

  return activeShift;
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::InitProblem()
{
  mQ.clear();
  mZoomRatios.clear();
  mShift.clear();
  mBoundaryShift.clear();
  mImprovementCoefficients.clear();

  if (mActiveConstraintNumber == 0)
    mActiveConstraintNumber = mConstraintIndeces.size();

  for (uint j = 0; j < mConstraintIndeces.size(); j++)
  {
    mImprovementCoefficients.push_back(10.0);
  }

  /// Сдвиг ограничений, оно же RHS
  mQ.resize(mConstraintIndeces.size());

  /// Коэфициенты масштабирования для функций ограничений
  mZoomRatios.resize(mConstraintIndeces.size());
  /// Покоордигатный сдвиг ограничений
  mShift.resize(mConstraintIndeces.size());

  for (uint i = 0; i < mConstraintIndeces.size(); i++)
  {
    mZoomRatios[i] = 1.0;
    mQ[i] = 0.0;
    mShift[i].resize(mDimension);
    for (int j = 0; j < mDimension; j++)
    {
      mShift[i][j] = 0.0;
    }
  }

  /// Покоордигатный сдвиг глобального минимума на границу
  mBoundaryShift.resize(mDimension);
  for (int i = 0; i < mDimension; i++)
  {
    mBoundaryShift[i] = 0.0;
  }

  /// Проверка параметров
  if (mProblemType != cptNormal)
  {
    for (uint j = 0; j < mFunctions.size(); j++)
    {
      if (!mFunctions[j].mIsMinimumKnown)
      {
        throw string("All global points must be known");
      }
    }
  }

  /// Задаем масштабирование для ограничений
  if (mIsZoom)
  {
    SetZoom();
  }

  /// Сдвигаем глобальные минимумы ограничений в координаты глобального минимума критерия
  if (mIsShift)
  {
    SetShift();
  }

  /// Изменяем ограничения для формирования заданной допустимой области
  if (mIsTotalDelta)
  {
    double q = CalculateRHS(mFeasibleDomainFraction);
    for (uint i = 0; i < mConstraintIndeces.size(); i++)
    {
      mQ[i] = q;
    }
  }

  /// Сдвигаем глобальный минимум критерия на границу допустимой области
  if (mIsBoundaryShift)
  {
    vector<double> shift(mConstraintIndeces.size());
    shift = SetBoundaryShift();
    for (uint i = 0; i < mConstraintIndeces.size(); i++)
    {
      mQ[i] = mQ[i] + shift[i];
    }
  }
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::TransformValue(double val, vector<double> point, int index) const
{
  double res = val;
  double resultCoefficient = 0;

  if (mIsImprovementOfTheObjective)
  {
    for (uint j = 0; j < mConstraintIndeces.size(); j++)
    {
      double val = Compute(mConstraintIndeces[j],
        TransformConstraintPoint(point, j));
      val = TransformConstraintValue(val, j);


      double fVal = std::max(val, 0.0);
      resultCoefficient += mImprovementCoefficients[j] * (fVal * fVal * fVal);
    }
  }

  return res - resultCoefficient;
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::TransformPoint(vector<double> point, int index) const
{
  vector<double> res(mDimension);

  if (mIsBoundaryShift)
  {
    double boundaryShift = 0; // величина сдвига
    vector<double> objectiveMin(mDimension);
    for (int j = 0; j < mDimension; j++)
    {
      objectiveMin[j] = mOptimumPoint[j];
      res[j] = point[j];
    }

    int coordinateNum = 0; // координата, по которой происходит сдвиг
    for (int i = 0; i < mDimension; i++)
    {
      if (mBoundaryShift[i] != 0)
      {
        coordinateNum = i;
        boundaryShift = mBoundaryShift[i];
      }
    }
    // преобразование координат таким образом, чтобы точка оптимума оказалась на границе
    // допустимой области
    if (boundaryShift > 0)
    {
      if ((point[coordinateNum] >= objectiveMin[coordinateNum]) &&
        (point[coordinateNum] < objectiveMin[coordinateNum] + boundaryShift))
      {
        res[coordinateNum] = objectiveMin[coordinateNum] + boundaryShift -
          (objectiveMin[coordinateNum] + boundaryShift - point[coordinateNum]) * 2;
      }
      else if ((point[coordinateNum] > objectiveMin[coordinateNum] - 2 * boundaryShift) &&
        (point[coordinateNum] < objectiveMin[coordinateNum]))
      {
        res[coordinateNum] = objectiveMin[coordinateNum] - 2 * boundaryShift +
          (point[coordinateNum] - (objectiveMin[coordinateNum] - 2 * boundaryShift)) / 2;
      }
    }
    else if (boundaryShift < 0)
    {
      if ((point[coordinateNum] > objectiveMin[coordinateNum]) &&
        (point[coordinateNum] < objectiveMin[coordinateNum] - 2 * boundaryShift))
      {
        res[coordinateNum] = objectiveMin[coordinateNum] - 2 * boundaryShift -
          (objectiveMin[coordinateNum] - 2 * boundaryShift - point[coordinateNum]) / 2;
      }
      else if ((point[coordinateNum] > objectiveMin[coordinateNum] + boundaryShift) &&
        (point[coordinateNum] < objectiveMin[coordinateNum]))
      {
        res[coordinateNum] = objectiveMin[coordinateNum] + boundaryShift +
          (point[coordinateNum] - (objectiveMin[coordinateNum] + boundaryShift)) * 2;
      }
    }
  }
  else
  {
    res = point;
  }

  return res;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::TransformConstraintValue(double val, int index) const
{
  double res;

  res = val - mQ[index];

  return res;
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::TransformConstraintPoint(vector<double> point,
  int index) const
{
  vector<double> res(mDimension);

  for (int i = 0; i < mDimension; i++)
  {
    res[i] = (point[i] - mShift[index][i]) /
      mZoomRatios[index];
  }

  return res;
}

// ------------------------------------------------------------------------------------------------
IConstrainedOptProblem::IConstrainedOptProblem(int dim, vector<double> loBound,
  vector<double> upBound, int probIndex,
  EConstrainedProblemType problemType, double fraction, int activeConstrNum)
  : IGeneralOptProblem(dim, loBound, upBound, probIndex)
{
  mFunctions.push_back({ mDimension, {}, 0, {}, 0,
    -1, true, false, false, false });
  mFunctionNumber = mFunctions.size();
  mCriterionIndeces.push_back(0);
  mProblemType = problemType;
  mFeasibleDomainFraction = fraction;
  mActiveConstraintNumber = activeConstrNum;

  mBoundarySearchPrecision = 20;

  if (problemType == cptInFeasibleDomain)
  {
    mIsShift = true;
    mIsZoom = true;
    mIsBoundaryShift = false;
    mIsImprovementOfTheObjective = false;
    mIsTotalDelta = true;
  }
  else if (problemType == cptOutFeasibleDomain)
  {
    mIsShift = true;
    mIsZoom = true;
    mIsBoundaryShift = false;
    mIsImprovementOfTheObjective = true;
    mIsTotalDelta = true;
  }
  else if (problemType == cptOnFeasibleBorder)
  {
    mIsShift = true;
    mIsZoom = true;
    mIsBoundaryShift = true;
    mIsImprovementOfTheObjective = false;
    mIsTotalDelta = true;
  }
  else
  {
    mIsShift = false;
    mIsZoom = false;
    mIsBoundaryShift = false;
    mIsImprovementOfTheObjective = false;
    mIsTotalDelta = false;
  }
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::GetOptimumPoint() const
{
  if (mProblemType == cptNormal)
  {
    if (mOptimumPoint.empty())
      throw string("Minimum point is unknown");
    for (uint j = 0; j < mConstraintIndeces.size(); j++)
      if (Compute(mConstraintIndeces[j], mOptimumPoint) > 0)
        throw string("Minimum point is unknown");
  }
  return TransformPoint(IGeneralOptProblem::GetOptimumPoint(), 0);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetOptimumValue() const
{
  if (mProblemType == cptNormal)
  {
    if (mOptimumPoint.empty())
      throw string("Minimum point is unknown");
    for (uint j = 0; j < mConstraintIndeces.size(); j++)
      if (Compute(mConstraintIndeces[j], mOptimumPoint) > 0)
        throw string("Minimum point is unknown");
  }
  return TransformValue(IGeneralOptProblem::GetOptimumValue(),
    IGeneralOptProblem::GetOptimumPoint(), 0);
}

// ------------------------------------------------------------------------------------------------
bool IConstrainedOptProblem::GetStatus(enum EOptFunctionParameter param) const
{
  return IGeneralOptProblem::GetStatus(mCriterionIndeces[0], param);
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::GetMaxPoint() const
{
  vector<double> point = IGeneralOptProblem::GetMaxPoint(mCriterionIndeces[0]);
  return TransformPoint(point, 0);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetMaxValue() const
{
  double val = IGeneralOptProblem::GetMaxValue(mCriterionIndeces[0]);
  return TransformValue(val, GetMaxPoint(), 0);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetLipschitzConstant() const
{
  double lip = IGeneralOptProblem::GetLipschitzConstant(0);
  return lip;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::ComputeFunction(const vector<double>& y) const
{
  double val = Compute(mCriterionIndeces[0],
    TransformPoint(y, 0));
  return TransformValue(val, y, 0);
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::ComputeFunctionDerivatives(const vector<double>& y) const
{
  vector<double> res = ComputeDerivatives(mCriterionIndeces[0], y);
  // TODO Реализовать преобразование, исходя из значения mProblemType и параметров задачи
  return res;
}

// ------------------------------------------------------------------------------------------------
int IConstrainedOptProblem::GetConstraintsNumber() const
{
  return mConstraintIndeces.size();
}

// ------------------------------------------------------------------------------------------------
bool IConstrainedOptProblem::GetConstraintStatus(int index, enum EOptFunctionParameter param) const
{
  return IGeneralOptProblem::GetStatus(mConstraintIndeces[index], param);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetConstraintLipschitzConstant(int index) const
{
  double val = IGeneralOptProblem::GetLipschitzConstant(index);
  return val;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::ComputeConstraint(int index, const vector<double>& y) const
{
  double val = Compute(mConstraintIndeces[index],
    TransformConstraintPoint(y, index));
  return TransformConstraintValue(val, index);
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::ComputeConstraints(const vector<double>& y,
  EConstraintComputationType t, int &index) const
{
  vector<double> res;
  switch (t)
  {
  case cctAllConstraints:
    for (uint i = 0; i < mConstraintIndeces.size(); i++)
    {
      res.push_back(ComputeConstraint(i, y));
    }
    index = mConstraintIndeces.size() - 1;
    break;
  case cctIndexScheme:
    double tmp;
    uint i;
    for (i = 0; i < mConstraintIndeces.size(); i++)
    {
      tmp = ComputeConstraint(i, y);
      if (tmp >= 0)
      {
        res.push_back(tmp);
        break;
      }
    }
    index = i;
    break;
  default:
    throw string("Unknown type of constraints computation");
  }
  return res;
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::ComputeConstraintDerivatives(int index,
  const vector<double>& y) const
{
  vector<double> res = ComputeDerivatives(mConstraintIndeces[index], y);
  // TODO Реализовать преобразование, исходя из значения mProblemType и параметров задачи
  return res;
}

// ------------------------------------------------------------------------------------------------
IConstrainedOptProblem::~IConstrainedOptProblem()
{
}