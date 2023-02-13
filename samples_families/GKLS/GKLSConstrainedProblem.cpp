#define _USE_MATH_DEFINES

#include "GKLSConstrainedProblem.hpp"

#include <math.h>
#include <algorithm>
#include <cassert>

// #include <fstream>

// ------------------------------------------------------------------------------------------------
TGKLSConstrainedProblem::TGKLSConstrainedProblem(EConstrainedProblemType problemType, double fraction,
  int activeConstrNum, int problemIndex, int dim, GKLSClass type, GKLSFuncionType functionType)
  : IConstrainedOptProblem(2, {}, {}, -1, problemType,
    fraction, activeConstrNum)
{
  options.resize(3);

  mDimension = dim;
  if ((mDimension <= 1) || (mDimension >= NUM_RND))
    throw "GKLS_DIM_ERROR"; /* problem dimension error */
  mProblemIndex = problemIndex;
  mLoBound.resize(mDimension);
  mUpBound.resize(mDimension);
  for (int i = 0; i < mDimension; i++) {
    mLoBound[i] = -1.0;
    mUpBound[i] = 1.0;
  }

  options[0].isArgSet = 0;
  options[0].mIsGeneratorMemoryAllocated = false;
  options[0].mIsDomainMemeoryAllocated = true;
  options[0].rnd_num = new double[NUM_RND];
  options[0].rand_condition = new double[KK];
  options[0].mFunctionType = functionType;
  options[0].index = problemIndex;

  options[0].SetFunctionClass(type, mDimension);
  SetOptions(0);

  mOptimumPoint.resize(mDimension);

  if (options[0].isArgSet == 1)
  {
    double *x = options[0].GKLS_minima.local_min[options[0].GKLS_glob.gm_index[0]];
    for (int i = 0; i < mDimension; i++)
      mOptimumPoint[i] = x[i];
  }

  mOptimumValue = Compute(0, mOptimumPoint);

  mFunctions[0].mDimension = mDimension;
  mFunctions[0].mMinimumPoint = mOptimumPoint;
  mFunctions[0].mMinimumValue = mOptimumValue;
  mFunctions[0].mIsMinimumKnown = true;
  mFunctions[0].mIsDerivativesKnown = true;

  options[1].index = problemIndex % 100 + 1;
  options[1].isArgSet = 0;
  options[1].mIsGeneratorMemoryAllocated = false;
  options[1].mIsDomainMemeoryAllocated = true;
  options[1].rnd_num = new double[NUM_RND];
  options[1].rand_condition = new double[KK];
  options[1].mFunctionType = functionType;
  options[1].SetFunctionClass(type, mDimension);
  SetOptions(1);

  vector<double> constraintOptPoint(mDimension);
  if (options[1].isArgSet == 1)
  {
    double *x = options[1].GKLS_minima.local_min[options[1].GKLS_glob.gm_index[0]];
    for (int i = 0; i < mDimension; i++)
      constraintOptPoint[i] = x[i];
  }
  mFunctions.push_back({ mDimension, constraintOptPoint,
    Compute(1,constraintOptPoint), {}, 0, -1, true, false, false, true });

  options[2].index = (problemIndex % 100 + 1) % 100 + 1;
  options[2].isArgSet = 0;
  options[2].mIsGeneratorMemoryAllocated = false;
  options[2].mIsDomainMemeoryAllocated = true;
  options[2].rnd_num = new double[NUM_RND];
  options[2].rand_condition = new double[KK];
  options[2].mFunctionType = functionType;
  options[2].SetFunctionClass(type, mDimension);
  SetOptions(2);

  if (options[2].isArgSet == 1)
  {
    double *x = options[2].GKLS_minima.local_min[options[2].GKLS_glob.gm_index[0]];
    for (int i = 0; i < mDimension; i++)
      constraintOptPoint[i] = x[i];
  }
  mFunctions.push_back({ mDimension, constraintOptPoint,
    Compute(2,constraintOptPoint),{}, 0, -1, true, false, false, true });

  mFunctionNumber = mFunctions.size();

  mCriterionIndeces.push_back(0);

  mConstraintIndeces.push_back(1);
  mConstraintIndeces.push_back(2);

  InitProblem();
}

// ------------------------------------------------------------------------------------------------
TGKLSConstrainedProblem::~TGKLSConstrainedProblem()
{
  for (unsigned i = 0; i < options.size(); i++)
  {
    if (options[i].mIsGeneratorMemoryAllocated)
      options[i].GKLS_free();

    delete[] options[i].rnd_num;
    delete[] options[i].rand_condition;
  }
}

// ------------------------------------------------------------------------------------------------
void TGKLSConstrainedProblem::SetOptions(int problemIndex)
{
  if (options[problemIndex].mIsGeneratorMemoryAllocated)
    options[problemIndex].GKLS_free();
  int err_code = options[problemIndex].GKLS_arg_generate(mProblemIndex, 0, mDimension, mLoBound, mUpBound);
  assert(err_code == GKLS_OK);
}

// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::Compute(int problemIndex, const vector<double>& y) const
{
  double value = 0.;

  switch (options[problemIndex].mFunctionType)
  {
  case TND:
    value = CalculateNDFunction(problemIndex, y.data());
    break;
  case TD:
    value = CalculateDFunction(problemIndex, y.data());
    break;
  case TD2:
    value = CalculateD2Function(problemIndex, y.data());
  }

  return value;
}

// ------------------------------------------------------------------------------------------------
vector<double> TGKLSConstrainedProblem::ComputeDerivatives(int problemIndex, const vector<double>& y) const
{
  vector<double> value(mDimension, 0);

  switch (options[problemIndex].mFunctionType)
  {
  case TND:
    throw "Not Derivatives";
    break;
  case TD:
    for (int i = 0; i < mDimension; i++)
      value[i] = CalculateDFunctionDeriv(problemIndex, i + 1, y.data());
    break;
  case TD2:
    for (int i = 0; i < mDimension; i++)
      value[i] = CalculateD2FunctionDeriv1(problemIndex, i + 1, y.data());
  }

  return value;
}



// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::CalculateNDFunction(int problemIndex, const double* x) const
{
  int i;
  unsigned index;

  double norm, scal, a, rho; /* working variables */

  if (!options[problemIndex].isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < options[problemIndex].GKLS_num_minima) &&
    (options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension) > options[problemIndex].GKLS_minima.rho[index]))
    index++;
  if (index == options[problemIndex].GKLS_num_minima)
  {
    norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], x, mDimension);
    /* Return the value of the paraboloid function */
    return (norm * norm + options[problemIndex].GKLS_minima.f[0]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (options[problemIndex].GKLS_norm(x, options[problemIndex].GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return options[problemIndex].GKLS_minima.f[index];

  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], options[problemIndex].GKLS_minima.local_min[index], mDimension);
  a = norm * norm + options[problemIndex].GKLS_minima.f[0] - options[problemIndex].GKLS_minima.f[index];
  rho = options[problemIndex].GKLS_minima.rho[index];
  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - options[problemIndex].GKLS_minima.local_min[index][i]) *
    (options[problemIndex].GKLS_minima.local_min[0][i] - options[problemIndex].GKLS_minima.local_min[index][i]);
  /* Return the value of the quadratic interpolation function */
  return ((1.0 - 2.0 / rho * scal / norm + a / rho / rho)*norm*norm +
    options[problemIndex].GKLS_minima.f[index]);
}

// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::CalculateDFunction(int problemIndex, const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, a, rho; /* working variables */

  if (!options[problemIndex].isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < options[problemIndex].GKLS_num_minima) &&
    (options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension) > options[problemIndex].GKLS_minima.rho[index]))
    index++;
  if (index == options[problemIndex].GKLS_num_minima)
  {
    norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], x, mDimension);
    /* Return the value of the paraboloid function */
    return (norm * norm + options[problemIndex].GKLS_minima.f[0]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (options[problemIndex].GKLS_norm(x, options[problemIndex].GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return options[problemIndex].GKLS_minima.f[index];

  // std::ofstream ofstr("data_point.txt");
  // for (int k = 0; k < 3; k++) {
  // for (i = 0; i < options[k].GKLS_num_minima; i++) {
  //   ofstr << "LOC_MIN" << k << "[" << 2 * i + 1 << "] = " << options[k].GKLS_minima.local_min[i][0] << std::endl;
  //   ofstr << "LOC_MIN" << k << "[" << 2 * i + 2 << "] = " << options[k].GKLS_minima.local_min[i][1] << std::endl;
  // }
  // ofstr << std::endl;
  // for (i = 0; i < options[k].GKLS_num_minima; i++) {
  //   ofstr << "RHO" << k << "[" << i + 1 << "] = " << options[k].GKLS_minima.rho[i] << std::endl;
  // }
  // ofstr << std::endl;
  // for (i = 0; i < options[k].GKLS_num_minima; i++) {
  //   ofstr << "F" << k << "[" << i + 1 << "] = " << options[k].GKLS_minima.f[i] << std::endl;
  // }
  // }
  // ofstr.close();

  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], options[problemIndex].GKLS_minima.local_min[index], mDimension);
  a = norm * norm + options[problemIndex].GKLS_minima.f[0] - options[problemIndex].GKLS_minima.f[index];
  rho = options[problemIndex].GKLS_minima.rho[index];
  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - options[problemIndex].GKLS_minima.local_min[index][i]) *
    (options[problemIndex].GKLS_minima.local_min[0][i] - options[problemIndex].GKLS_minima.local_min[index][i]);
  /* Return the value of the cubic interpolation function */
  return (2.0 / rho / rho * scal / norm - 2.0*a / rho / rho / rho)*norm*norm*norm +
    (1.0 - 4.0*scal / norm / rho + 3.0*a / rho / rho)*norm*norm + options[problemIndex].GKLS_minima.f[index];
}

// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::CalculateD2Function(int problemIndex, const double* x) const
{
  unsigned int dim, i, index;
  double norm, scal, a, rho; /* working variables */

  if (!options[problemIndex].isArgSet) return GKLS_MAX_VALUE;

  dim = (unsigned)mDimension;
  for (i = 0; i < dim; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) ||
      (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima */
  index = 1;
  while ((index < options[problemIndex].GKLS_num_minima) &&
    (options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension) > options[problemIndex].GKLS_minima.rho[index]))
    index++;
  if (index == options[problemIndex].GKLS_num_minima)
  {
    norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], x, mDimension);
    /* Return the value of the paraboloid function */
    return (norm * norm + options[problemIndex].GKLS_minima.f[0]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (options[problemIndex].GKLS_norm(x, options[problemIndex].GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return options[problemIndex].GKLS_minima.f[index];

  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], options[problemIndex].GKLS_minima.local_min[index], mDimension);
  a = norm * norm + options[problemIndex].GKLS_minima.f[0] - options[problemIndex].GKLS_minima.f[index];
  rho = options[problemIndex].GKLS_minima.rho[index];
  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < dim; i++)
    scal += (x[i] - options[problemIndex].GKLS_minima.local_min[index][i]) *
    (options[problemIndex].GKLS_minima.local_min[0][i] - options[problemIndex].GKLS_minima.local_min[index][i]);
  /* Return the value of the quintic interpolation function */
  return ((-6.0*scal / norm / rho + 6.0*a / rho / rho + 1.0 - options[problemIndex].delta / 2.0) *
    norm * norm / rho / rho +
    (16.0*scal / norm / rho - 15.0*a / rho / rho - 3.0 + 1.5*options[problemIndex].delta) * norm / rho +
    (-12.0*scal / norm / rho + 10.0*a / rho / rho + 3.0 - 1.5*options[problemIndex].delta)) *
    norm * norm * norm / rho +
    0.5*options[problemIndex].delta*norm*norm + options[problemIndex].GKLS_minima.f[index];
}

// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::CalculateDFunctionDeriv(int problemIndex, unsigned var_j, const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, dif, a, rho, h; /* working variables */

  if ((var_j == 0) || (var_j > (unsigned)mDimension)) return GKLS_MAX_VALUE;
  else  var_j = var_j - 1; /* to be used as an index of array */

  if (!options[problemIndex].isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < options[problemIndex].GKLS_num_minima) &&
    (options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension) > options[problemIndex].GKLS_minima.rho[index]))
    index++;
  if (index == options[problemIndex].GKLS_num_minima)
  {
    /* Return the value of the first order partial derivative of the paraboloid function */
    return 2.0*(x[var_j] - options[problemIndex].GKLS_minima.local_min[0][var_j]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (options[problemIndex].GKLS_norm(x, options[problemIndex].GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return 0.0;

  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], options[problemIndex].GKLS_minima.local_min[index], mDimension);
  a = norm * norm + options[problemIndex].GKLS_minima.f[0] - options[problemIndex].GKLS_minima.f[index];
  rho = options[problemIndex].GKLS_minima.rho[index];
  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - options[problemIndex].GKLS_minima.local_min[index][i]) *
    (options[problemIndex].GKLS_minima.local_min[0][i] - options[problemIndex].GKLS_minima.local_min[index][i]);
  dif = x[var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j];
  h = (options[problemIndex].GKLS_minima.local_min[0][var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j])*norm -
    scal*dif / norm;
  /* Return the value of dC(x)/dx[var_i] of the D-type function */
  return (h * (2.0 / rho / rho*norm - 4.0 / rho) +
    dif * (6.0 / rho / rho*scal - 6.0 / rho / rho / rho*a*norm -
      8.0 / rho / norm*scal + 6.0 / rho / rho*a + 2.0));
}

// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::CalculateD2FunctionDeriv1(int problemIndex, unsigned var_j, const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, dif, a, rho, h; /* working variables */

  if ((var_j == 0) || (var_j > (unsigned)mDimension)) return GKLS_MAX_VALUE;
  else  var_j = var_j - 1; /* to be used as an index of array */

  if (!options[problemIndex].isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < options[problemIndex].GKLS_num_minima) &&
    (options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension) > options[problemIndex].GKLS_minima.rho[index]))
    index++;
  if (index == options[problemIndex].GKLS_num_minima)
  {
    /* Return the value of the first order partial derivative of the paraboloid function */
    return 2.0*(x[var_j] - options[problemIndex].GKLS_minima.local_min[0][var_j]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (options[problemIndex].GKLS_norm(x, options[problemIndex].GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return 0.0;

  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], options[problemIndex].GKLS_minima.local_min[index], mDimension);
  a = norm * norm + options[problemIndex].GKLS_minima.f[0] - options[problemIndex].GKLS_minima.f[index];
  rho = options[problemIndex].GKLS_minima.rho[index];
  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - options[problemIndex].GKLS_minima.local_min[index][i]) *
    (options[problemIndex].GKLS_minima.local_min[0][i] - options[problemIndex].GKLS_minima.local_min[index][i]);
  dif = x[var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j];
  h = (options[problemIndex].GKLS_minima.local_min[0][var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j])*norm -
    scal*dif / norm;
  /* Return the value of dQ(x)/dx[var_i] of the D2-type function */
  return (h*norm / rho / rho * (-6.0*norm*norm / rho / rho + 16.0*norm / rho - 12.0) +
    dif*norm * ((-30.0 / rho / norm*scal + 30 / rho / rho*a + 5.0 - 2.5*options[problemIndex].delta) / rho / rho / rho*norm*norm +
    (64.0 / rho / norm*scal - 60.0 / rho / rho*a - 12.0 + 6.0*options[problemIndex].delta) / rho / rho*norm +
      (-36.0 / rho / norm*scal + 30.0 / rho / rho*a + 9.0 - 4.5*options[problemIndex].delta) / rho) +
    dif*options[problemIndex].delta);
}

// ------------------------------------------------------------------------------------------------
double TGKLSConstrainedProblem::CalculateD2FunctionDeriv2(int problemIndex, unsigned var_j, unsigned var_k, const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, a, rho,
    dh, difj, difk, hj, hk, dQ_jk; /* working variables */
  int the_same;  /* is TRUE if var_j==var_k */

  if ((var_j == 0) || (var_j > (unsigned)mDimension)) return GKLS_MAX_VALUE;
  if ((var_k == 0) || (var_k > (unsigned)mDimension)) return GKLS_MAX_VALUE;
  the_same = (var_j == var_k);
  var_j = var_j - 1; var_k = var_k - 1; /* to be used as indexes of array */

  if (!options[problemIndex].isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < options[problemIndex].GKLS_num_minima) &&
    (options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension) > options[problemIndex].GKLS_minima.rho[index]))
    index++;
  if (index == options[problemIndex].GKLS_num_minima)
  {
    /* Return the value of the second order partial derivative of the paraboloid function */
    if (the_same) return 2.0;
    else return 0.0;
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (options[problemIndex].GKLS_norm(x, options[problemIndex].GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
  {
    if (the_same) return options[problemIndex].delta;
    else return 0.0;
  }

  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[0], options[problemIndex].GKLS_minima.local_min[index], mDimension);
  a = norm * norm + options[problemIndex].GKLS_minima.f[0] - options[problemIndex].GKLS_minima.f[index];
  rho = options[problemIndex].GKLS_minima.rho[index];
  norm = options[problemIndex].GKLS_norm(options[problemIndex].GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - options[problemIndex].GKLS_minima.local_min[index][i]) *
    (options[problemIndex].GKLS_minima.local_min[0][i] - options[problemIndex].GKLS_minima.local_min[index][i]);
  difj = x[var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j];
  difk = x[var_k] - options[problemIndex].GKLS_minima.local_min[index][var_k];
  hj = (options[problemIndex].GKLS_minima.local_min[0][var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j])*norm -
    scal*difj / norm;
  hk = (options[problemIndex].GKLS_minima.local_min[0][var_k] - options[problemIndex].GKLS_minima.local_min[index][var_k])*norm -
    scal*difk / norm;

  dh = (options[problemIndex].GKLS_minima.local_min[0][var_j] - options[problemIndex].GKLS_minima.local_min[index][var_j])*difk / norm -
    hk*difj / norm / norm;
  if (the_same) dh = dh - scal / norm;

  dQ_jk = -6.0 / rho / rho / rho / rho*(dh*norm*norm*norm + 3.0*hj*difk*norm) -
    30.0 / rho / rho / rho / rho*hk*difj*norm +
    15.0 / rho / rho / rho*(-6.0 / rho*scal / norm + 6.0 / rho / rho*a + 1 - 0.5*options[problemIndex].delta)*difj*difk*norm +
    16.0 / rho / rho / rho*(dh*norm*norm + 2.0*hj*difk) +
    64.0 / rho / rho / rho*hk*difj +
    8.0 / rho / rho*(16.0 / rho*scal / norm - 15.0 / rho / rho*a - 3.0 + 1.5*options[problemIndex].delta)*difj*difk -
    12.0 / rho / rho*(dh*norm + hj*difk / norm) -
    36.0 / rho / rho*hk*difj / norm +
    3.0 / rho*(-12.0 / rho*scal / norm + 10.0 / rho / rho*a + 3.0 - 1.5*options[problemIndex].delta)*difj*difk / norm;

  if (the_same)
    dQ_jk = dQ_jk +
    5.0*norm*norm*norm / rho / rho / rho*(-6.0 / rho*scal / norm + 6.0 / rho / rho*a + 1 - 0.5*options[problemIndex].delta) +
    4.0*norm*norm / rho / rho*(16.0 / rho*scal / norm - 15.0 / rho / rho*a - 3.0 + 1.5*options[problemIndex].delta) +
    3.0*norm / rho*(-12.0 / rho*scal / norm + 10.0 / rho / rho*a + 3.0 - 1.5*options[problemIndex].delta) +
    options[problemIndex].delta;
  /* Return the value of d^2[Q(x)]/dx[var_j]dx[var_k] of the D2-type function */
  return dQ_jk;
}

// ------------------------------------------------------------------------------------------------
int TGKLSConstrainedProblem::CalculateDFunctionGradient(int problemIndex, const double* x, double* g) const
{
  int i;
  int error_code = GKLS_OK;

  if (!options[problemIndex].isArgSet) return GKLS_DERIV_EVAL_ERROR;
  if (g == NULL) return GKLS_DERIV_EVAL_ERROR;

  for (i = 1; i <= mDimension; i++)
  {
    g[i - 1] = CalculateDFunctionDeriv(problemIndex, i, x);
    if (g[i - 1] >= GKLS_MAX_VALUE - 1000.0)
      error_code = GKLS_DERIV_EVAL_ERROR;
  }
  return error_code;
}

// ------------------------------------------------------------------------------------------------
int TGKLSConstrainedProblem::CalculateD2FunctionGradient(int problemIndex, const double* x, double* g) const
{
  int i;
  int error_code = GKLS_OK;

  if (!options[problemIndex].isArgSet) return GKLS_DERIV_EVAL_ERROR;
  if (g == NULL) return GKLS_DERIV_EVAL_ERROR;

  for (i = 1; i <= mDimension; i++)
  {
    g[i - 1] = CalculateD2FunctionDeriv1(problemIndex, i, x);
    if (g[i - 1] >= GKLS_MAX_VALUE - 1000.0)
      error_code = GKLS_DERIV_EVAL_ERROR;
  }
  return error_code;
}

// ------------------------------------------------------------------------------------------------
int TGKLSConstrainedProblem::CalculateD2FunctionHessian(int problemIndex, const double* x, double** h) const
{
  int i, j;
  int error_code = GKLS_OK;

  if (!options[problemIndex].isArgSet) return GKLS_DERIV_EVAL_ERROR;
  if (h == NULL) return GKLS_DERIV_EVAL_ERROR;
  for (i = 1; i <= mDimension; i++)
    if (h[i - 1] == NULL) return GKLS_DERIV_EVAL_ERROR;

  for (i = 1; i <= mDimension; i++)
    for (j = 1; j <= mDimension; j++)
    {
      h[i - 1][j - 1] = CalculateD2FunctionDeriv2(problemIndex, i, j, x);
      if (h[i - 1][j - 1] >= GKLS_MAX_VALUE - 1000.0)
        error_code = GKLS_DERIV_EVAL_ERROR;
    }
  return error_code;
}