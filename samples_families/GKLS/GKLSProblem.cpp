#define _USE_MATH_DEFINES

#include <math.h>
#include <algorithm>
#include <cassert>

#include "GKLSProblem.hpp"

// #include <fstream>

// ------------------------------------------------------------------------------------------------
TGKLSProblem::TGKLSProblem(int problemIndex, int dim, GKLSClass type, GKLSFuncionType functionType) : IOptProblem()
{
  mDimension = dim;
  mProblemIndex = problemIndex;
  mLoBound.resize(mDimension);
  mUpBound.resize(mDimension);
  for (int i = 0; i < mDimension; i++) {
    mLoBound[i] = -1.0;
    mUpBound[i] = 1.0;
  }

  option.isArgSet = 0;
  option.mIsGeneratorMemoryAllocated = false;
  option.mIsDomainMemeoryAllocated = true;
  option.rnd_num = new double[NUM_RND];
  option.rand_condition = new double[KK];
  option.mFunctionType = functionType;
  option.index = problemIndex;

  option.SetFunctionClass(type, mDimension);
  SetOptions();

  if ((mDimension <= 1) || (mDimension >= NUM_RND))
    throw "GKLS_DIM_ERROR"; /* problem dimension error */


  mOptimumPoint.resize(mDimension);

  if (option.isArgSet == 1)
  {
    double *x = option.GKLS_minima.local_min[option.GKLS_glob.gm_index[0]];
    for (int i = 0; i < mDimension; i++)
      mOptimumPoint[i] = x[i];
  }

  mOptimumValue = Compute(0, mOptimumPoint);
  TOptFunction func;
  func.mDimension = mDimension;
  func.mMinimumPoint = mOptimumPoint;
  func.mMinimumValue = mOptimumValue;
  func.mIsMinimumKnown = true;
  func.mIsDerivativesKnown = true;
  mFunctions.push_back(func);
  mFunctionNumber = 1;
  mFunctionIndex = 0;
}

// ------------------------------------------------------------------------------------------------
TGKLSProblem::~TGKLSProblem()
{
  if (option.mIsGeneratorMemoryAllocated)
    option.GKLS_free();

  delete[] option.rnd_num;
  delete[] option.rand_condition;
}

// ------------------------------------------------------------------------------------------------
void TGKLSProblem::SetOptions()
{
  if (option.mIsGeneratorMemoryAllocated)
    option.GKLS_free();
  int err_code = option.GKLS_arg_generate(mProblemIndex, 0, mDimension, mLoBound, mUpBound);
  assert(err_code == GKLS_OK);

}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::Compute(int index, const vector<double>& y) const
{
  double value = 0.;

  switch (option.mFunctionType)
  {
  case TND:
    value = CalculateNDFunction(y.data());
    break;
  case TD:
    value = CalculateDFunction(y.data());
    break;
  case TD2:
    value = CalculateD2Function(y.data());
  }

  return value;
}

// ------------------------------------------------------------------------------------------------
vector<double> TGKLSProblem::ComputeDerivatives(int index, const vector<double>& y) const
{
  vector<double> value(mDimension, 0);

  switch (option.mFunctionType)
  {
  case TND:
    throw "Not Derivatives";
    break;
  case TD:
    for (int i = 0; i < mDimension; i++)
      value[i] = CalculateDFunctionDeriv(i + 1, y.data());
    break;
  case TD2:
    for (int i = 0; i < mDimension; i++)
      value[i] = CalculateD2FunctionDeriv1(i + 1, y.data());
  }

  return value;
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateNDFunction(const double* x) const
{
  int i;
  unsigned index;

  double norm, scal, a, rho; /* working variables */

  if (!option.isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < option.GKLS_num_minima) &&
    (option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension) > option.GKLS_minima.rho[index]))
    index++;
  if (index == option.GKLS_num_minima)
  {
    norm = option.GKLS_norm(option.GKLS_minima.local_min[0], x, mDimension);
    /* Return the value of the paraboloid function */
    return (norm * norm + option.GKLS_minima.f[0]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (option.GKLS_norm(x, option.GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return option.GKLS_minima.f[index];

  norm = option.GKLS_norm(option.GKLS_minima.local_min[0], option.GKLS_minima.local_min[index], mDimension);
  a = norm * norm + option.GKLS_minima.f[0] - option.GKLS_minima.f[index];
  rho = option.GKLS_minima.rho[index];
  norm = option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - option.GKLS_minima.local_min[index][i]) *
    (option.GKLS_minima.local_min[0][i] - option.GKLS_minima.local_min[index][i]);
  /* Return the value of the quadratic interpolation function */
  return ((1.0 - 2.0 / rho * scal / norm + a / rho / rho)*norm*norm +
    option.GKLS_minima.f[index]);
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateDFunction(const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, a, rho; /* working variables */

  if (!option.isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */

  // std::ofstream ofstr("data_point.txt");
  // for (i = 0; i < option.GKLS_num_minima; i++) {
  //   ofstr << "LOC_MIN[" << 2 * i + 1 << "] = " << option.GKLS_minima.local_min[i][0] << std::endl;
  //   ofstr << "LOC_MIN[" << 2 * i + 2 << "] = " << option.GKLS_minima.local_min[i][1] << std::endl;
  // }
  // ofstr << std::endl;
  // for (i = 0; i < option.GKLS_num_minima; i++) {
  //   ofstr << "RHO[" << i + 1 << "] = " << option.GKLS_minima.rho[i] << std::endl;
  // }
  // ofstr << std::endl;
  // for (i = 0; i < option.GKLS_num_minima; i++) {
  //   ofstr << "F[" << i + 1 << "] = " << option.GKLS_minima.f[i] << std::endl;
  // }
  // ofstr.close();

  index = 1;
  while ((index < option.GKLS_num_minima) &&
    (option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension) > option.GKLS_minima.rho[index]))
    index++;
  if (index == option.GKLS_num_minima)
  {
    norm = option.GKLS_norm(option.GKLS_minima.local_min[0], x, mDimension);
    /* Return the value of the paraboloid function */
    return (norm * norm + option.GKLS_minima.f[0]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (option.GKLS_norm(x, option.GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return option.GKLS_minima.f[index];

  norm = option.GKLS_norm(option.GKLS_minima.local_min[0], option.GKLS_minima.local_min[index], mDimension);
  a = norm * norm + option.GKLS_minima.f[0] - option.GKLS_minima.f[index];
  rho = option.GKLS_minima.rho[index];
  norm = option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - option.GKLS_minima.local_min[index][i]) *
    (option.GKLS_minima.local_min[0][i] - option.GKLS_minima.local_min[index][i]);
  /* Return the value of the cubic interpolation function */
  return (2.0 / rho / rho * scal / norm - 2.0*a / rho / rho / rho)*norm*norm*norm +
    (1.0 - 4.0*scal / norm / rho + 3.0*a / rho / rho)*norm*norm + option.GKLS_minima.f[index];
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateD2Function(const double* x) const
{
  unsigned int dim, i, index;
  double norm, scal, a, rho; /* working variables */

  if (!option.isArgSet) return GKLS_MAX_VALUE;

  dim = (unsigned)mDimension;
  for (i = 0; i < dim; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) ||
      (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima */
  index = 1;
  while ((index < option.GKLS_num_minima) &&
    (option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension) > option.GKLS_minima.rho[index]))
    index++;
  if (index == option.GKLS_num_minima)
  {
    norm = option.GKLS_norm(option.GKLS_minima.local_min[0], x, mDimension);
    /* Return the value of the paraboloid function */
    return (norm * norm + option.GKLS_minima.f[0]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (option.GKLS_norm(x, option.GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return option.GKLS_minima.f[index];

  norm = option.GKLS_norm(option.GKLS_minima.local_min[0], option.GKLS_minima.local_min[index], mDimension);
  a = norm * norm + option.GKLS_minima.f[0] - option.GKLS_minima.f[index];
  rho = option.GKLS_minima.rho[index];
  norm = option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < dim; i++)
    scal += (x[i] - option.GKLS_minima.local_min[index][i]) *
    (option.GKLS_minima.local_min[0][i] - option.GKLS_minima.local_min[index][i]);
  /* Return the value of the quintic interpolation function */
  return ((-6.0*scal / norm / rho + 6.0*a / rho / rho + 1.0 - option.delta / 2.0) *
    norm * norm / rho / rho +
    (16.0*scal / norm / rho - 15.0*a / rho / rho - 3.0 + 1.5*option.delta) * norm / rho +
    (-12.0*scal / norm / rho + 10.0*a / rho / rho + 3.0 - 1.5*option.delta)) *
    norm * norm * norm / rho +
    0.5*option.delta*norm*norm + option.GKLS_minima.f[index];
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateDFunctionDeriv(unsigned var_j, const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, dif, a, rho, h; /* working variables */

  if ((var_j == 0) || (var_j > (unsigned)mDimension)) return GKLS_MAX_VALUE;
  else  var_j = var_j - 1; /* to be used as an index of array */

  if (!option.isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < option.GKLS_num_minima) &&
    (option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension) > option.GKLS_minima.rho[index]))
    index++;
  if (index == option.GKLS_num_minima)
  {
    /* Return the value of the first order partial derivative of the paraboloid function */
    return 2.0*(x[var_j] - option.GKLS_minima.local_min[0][var_j]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (option.GKLS_norm(x, option.GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return 0.0;

  norm = option.GKLS_norm(option.GKLS_minima.local_min[0], option.GKLS_minima.local_min[index], mDimension);
  a = norm * norm + option.GKLS_minima.f[0] - option.GKLS_minima.f[index];
  rho = option.GKLS_minima.rho[index];
  norm = option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - option.GKLS_minima.local_min[index][i]) *
    (option.GKLS_minima.local_min[0][i] - option.GKLS_minima.local_min[index][i]);
  dif = x[var_j] - option.GKLS_minima.local_min[index][var_j];
  h = (option.GKLS_minima.local_min[0][var_j] - option.GKLS_minima.local_min[index][var_j])*norm -
    scal*dif / norm;
  /* Return the value of dC(x)/dx[var_i] of the D-type function */
  return (h * (2.0 / rho / rho*norm - 4.0 / rho) +
    dif * (6.0 / rho / rho*scal - 6.0 / rho / rho / rho*a*norm -
      8.0 / rho / norm*scal + 6.0 / rho / rho*a + 2.0));
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateD2FunctionDeriv1(unsigned var_j, const double* x) const
{
  int i;
  unsigned index;
  double norm, scal, dif, a, rho, h; /* working variables */

  if ((var_j == 0) || (var_j > (unsigned) mDimension)) return GKLS_MAX_VALUE;
  else  var_j = var_j - 1; /* to be used as an index of array */

  if (!option.isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < option.GKLS_num_minima) &&
    (option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension) > option.GKLS_minima.rho[index]))
    index++;
  if (index == option.GKLS_num_minima)
  {
    /* Return the value of the first order partial derivative of the paraboloid function */
    return 2.0*(x[var_j] - option.GKLS_minima.local_min[0][var_j]);
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (option.GKLS_norm(x, option.GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
    return 0.0;

  norm = option.GKLS_norm(option.GKLS_minima.local_min[0], option.GKLS_minima.local_min[index], mDimension);
  a = norm * norm + option.GKLS_minima.f[0] - option.GKLS_minima.f[index];
  rho = option.GKLS_minima.rho[index];
  norm = option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - option.GKLS_minima.local_min[index][i]) *
    (option.GKLS_minima.local_min[0][i] - option.GKLS_minima.local_min[index][i]);
  dif = x[var_j] - option.GKLS_minima.local_min[index][var_j];
  h = (option.GKLS_minima.local_min[0][var_j] - option.GKLS_minima.local_min[index][var_j])*norm -
    scal*dif / norm;
  /* Return the value of dQ(x)/dx[var_i] of the D2-type function */
  return (h*norm / rho / rho * (-6.0*norm*norm / rho / rho + 16.0*norm / rho - 12.0) +
    dif*norm * ((-30.0 / rho / norm*scal + 30 / rho / rho*a + 5.0 - 2.5*option.delta) / rho / rho / rho*norm*norm +
    (64.0 / rho / norm*scal - 60.0 / rho / rho*a - 12.0 + 6.0*option.delta) / rho / rho*norm +
      (-36.0 / rho / norm*scal + 30.0 / rho / rho*a + 9.0 - 4.5*option.delta) / rho) +
    dif*option.delta);
}

// ------------------------------------------------------------------------------------------------
double TGKLSProblem::CalculateD2FunctionDeriv2(unsigned var_j, unsigned var_k, const double* x) const
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

  if (!option.isArgSet) return GKLS_MAX_VALUE;

  for (i = 0; i < mDimension; i++)
    if ((x[i] < mLoBound[i] - GKLS_PRECISION) || (x[i] > mUpBound[i] + GKLS_PRECISION))
      return GKLS_MAX_VALUE;
  /* Check wether x belongs to some basin of local minima, M(index) <> T */
  /* Attention: number of local minima must be >= 2 */
  index = 1;
  while ((index < option.GKLS_num_minima) &&
    (option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension) > option.GKLS_minima.rho[index]))
    index++;
  if (index == option.GKLS_num_minima)
  {
    /* Return the value of the second order partial derivative of the paraboloid function */
    if (the_same) return 2.0;
    else return 0.0;
  }

  /* Check wether x coincides with the local minimizer M(index) */
  if (option.GKLS_norm(x, option.GKLS_minima.local_min[index], mDimension) < GKLS_PRECISION)
  {
    if (the_same) return option.delta;
    else return 0.0;
  }

  norm = option.GKLS_norm(option.GKLS_minima.local_min[0], option.GKLS_minima.local_min[index], mDimension);
  a = norm * norm + option.GKLS_minima.f[0] - option.GKLS_minima.f[index];
  rho = option.GKLS_minima.rho[index];
  norm = option.GKLS_norm(option.GKLS_minima.local_min[index], x, mDimension);
  scal = 0.0;
  for (i = 0; i < mDimension; i++)
    scal += (x[i] - option.GKLS_minima.local_min[index][i]) *
    (option.GKLS_minima.local_min[0][i] - option.GKLS_minima.local_min[index][i]);
  difj = x[var_j] - option.GKLS_minima.local_min[index][var_j];
  difk = x[var_k] - option.GKLS_minima.local_min[index][var_k];
  hj = (option.GKLS_minima.local_min[0][var_j] - option.GKLS_minima.local_min[index][var_j])*norm -
    scal*difj / norm;
  hk = (option.GKLS_minima.local_min[0][var_k] - option.GKLS_minima.local_min[index][var_k])*norm -
    scal*difk / norm;

  dh = (option.GKLS_minima.local_min[0][var_j] - option.GKLS_minima.local_min[index][var_j])*difk / norm -
    hk*difj / norm / norm;
  if (the_same) dh = dh - scal / norm;

  dQ_jk = -6.0 / rho / rho / rho / rho*(dh*norm*norm*norm + 3.0*hj*difk*norm) -
    30.0 / rho / rho / rho / rho*hk*difj*norm +
    15.0 / rho / rho / rho*(-6.0 / rho*scal / norm + 6.0 / rho / rho*a + 1 - 0.5*option.delta)*difj*difk*norm +
    16.0 / rho / rho / rho*(dh*norm*norm + 2.0*hj*difk) +
    64.0 / rho / rho / rho*hk*difj +
    8.0 / rho / rho*(16.0 / rho*scal / norm - 15.0 / rho / rho*a - 3.0 + 1.5*option.delta)*difj*difk -
    12.0 / rho / rho*(dh*norm + hj*difk / norm) -
    36.0 / rho / rho*hk*difj / norm +
    3.0 / rho*(-12.0 / rho*scal / norm + 10.0 / rho / rho*a + 3.0 - 1.5*option.delta)*difj*difk / norm;

  if (the_same)
    dQ_jk = dQ_jk +
    5.0*norm*norm*norm / rho / rho / rho*(-6.0 / rho*scal / norm + 6.0 / rho / rho*a + 1 - 0.5*option.delta) +
    4.0*norm*norm / rho / rho*(16.0 / rho*scal / norm - 15.0 / rho / rho*a - 3.0 + 1.5*option.delta) +
    3.0*norm / rho*(-12.0 / rho*scal / norm + 10.0 / rho / rho*a + 3.0 - 1.5*option.delta) +
    option.delta;
  /* Return the value of d^2[Q(x)]/dx[var_j]dx[var_k] of the D2-type function */
  return dQ_jk;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::CalculateDFunctionGradient(const double* x, double* g) const
{
  int i;
  int error_code = GKLS_OK;

  if (!option.isArgSet) return GKLS_DERIV_EVAL_ERROR;
  if (g == NULL) return GKLS_DERIV_EVAL_ERROR;

  for (i = 1; i <= mDimension; i++)
  {
    g[i - 1] = CalculateDFunctionDeriv(i, x);
    if (g[i - 1] >= GKLS_MAX_VALUE - 1000.0)
      error_code = GKLS_DERIV_EVAL_ERROR;
  }
  return error_code;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::CalculateD2FunctionGradient(const double* x, double* g) const
{
  int i;
  int error_code = GKLS_OK;

  if (!option.isArgSet) return GKLS_DERIV_EVAL_ERROR;
  if (g == NULL) return GKLS_DERIV_EVAL_ERROR;

  for (i = 1; i <= mDimension; i++)
  {
    g[i - 1] = CalculateD2FunctionDeriv1(i, x);
    if (g[i - 1] >= GKLS_MAX_VALUE - 1000.0)
      error_code = GKLS_DERIV_EVAL_ERROR;
  }
  return error_code;
}

// ------------------------------------------------------------------------------------------------
int TGKLSProblem::CalculateD2FunctionHessian(const double* x, double** h) const
{
  int i, j;
  int error_code = GKLS_OK;

  if (!option.isArgSet) return GKLS_DERIV_EVAL_ERROR;
  if (h == NULL) return GKLS_DERIV_EVAL_ERROR;
  for (i = 1; i <= mDimension; i++)
    if (h[i - 1] == NULL) return GKLS_DERIV_EVAL_ERROR;

  for (i = 1; i <= mDimension; i++)
    for (j = 1; j <= mDimension; j++)
    {
      h[i - 1][j - 1] = CalculateD2FunctionDeriv2(i, j, x);
      if (h[i - 1][j - 1] >= GKLS_MAX_VALUE - 1000.0)
        error_code = GKLS_DERIV_EVAL_ERROR;
    }
  return error_code;
}