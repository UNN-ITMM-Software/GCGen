#define _USE_MATH_DEFINES

#include "GrishaginConstrainedProblem.hpp"

#include <math.h>
#include <algorithm>
#include <cmath>
#include <cassert>

// #include <fstream>

// ------------------------------------------------------------------------------------------------
GrishaginConstrainedProblem::GrishaginConstrainedProblem(EConstrainedProblemType problemType, double fraction,
  int activeConstrNum, int problemIndex)
  : IConstrainedOptProblem(2, {}, {}, -1, problemType,
    fraction, activeConstrNum)
{
  options.resize(3);
  mProblemIndex = problemIndex;
  options[0].index = problemIndex;
  SetOptions(0);

  mDimension = 2;
  mLoBound = { 0.0, 0.0 };
  mUpBound = { 1.0, 1.0 };

  mOptimumPoint = { rand_minimums[2 * (mProblemIndex - 1)], rand_minimums[2 * (mProblemIndex - 1) + 1] };
  mOptimumValue = Compute(0, mOptimumPoint);

  mFunctions[0].mDimension = mDimension;
  mFunctions[0].mMinimumPoint = mOptimumPoint;
  mFunctions[0].mMinimumValue = mOptimumValue;
  mFunctions[0].mIsMinimumKnown = true;
  mFunctions[0].mIsDerivativesKnown = true;

  options[1].index = problemIndex % 100 + 1;
  SetOptions(1);
  mFunctions.push_back({ mDimension,{ rand_minimums[2 * (options[1].index - 1)],
    rand_minimums[2 * (options[1].index - 1) + 1] },
    Compute(1,{ rand_minimums[2 * (options[1].index - 1)],
      rand_minimums[2 * (options[1].index - 1) + 1] }),
      {}, 0, -1,
    true, false, false, true });

  options[2].index = (problemIndex % 100 + 1) % 100 + 1;
  SetOptions(2);
  mFunctions.push_back({ mDimension,{ rand_minimums[2 * (options[2].index - 1)],
    rand_minimums[2 * (options[2].index - 1) + 1] },
    Compute(1,{ rand_minimums[2 * (options[2].index - 1)],
      rand_minimums[2 * (options[2].index - 1) + 1] }),
      {}, 0, -1,
    true, false, false, true });

  mFunctionNumber = mFunctions.size();

  mCriterionIndeces.push_back(0);

  mConstraintIndeces.push_back(1);
  mConstraintIndeces.push_back(2);

  InitProblem();
}

// ------------------------------------------------------------------------------------------------
void GrishaginConstrainedProblem::SetOptions(int index)
{
  assert(mProblemIndex > 0 && mProblemIndex <= 100);

  int lst, i, j, i1, i2, i3;
  int nf = options[index].index;

  if (nf < 1 || nf>100)
    nf = 1;
  lst = 10;
  i1 = (nf - 1) / lst;
  i2 = i1*lst;
  for (j = 0; j < 45; j++)
    options[index].icnf[j] = matcon[i1][j];
  if (i2 != (nf - 1)) {
    i3 = nf - 1 - i2;
    for (j = 1; j <= i3; j++)
      for (i = 0; i < 196; i++)
        rndm20(options[index].icnf);
  }
  for (j = 0; j < 7; j++)
    for (i = 0; i < 7; i++) {
      options[index].af[i][j] = 2.*rndm20(options[index].icnf) - 1.;
      options[index].cf[i][j] = 2.*rndm20(options[index].icnf) - 1.;
    }
  for (j = 0; j < 7; j++)
    for (i = 0; i < 7; i++) {
      options[index].bf[i][j] = 2.*rndm20(options[index].icnf) - 1.;
      options[index].df[i][j] = 2.*rndm20(options[index].icnf) - 1.;
    }
}

// ------------------------------------------------------------------------------------------------
double GrishaginConstrainedProblem::Compute(int index, const vector<double>& y) const
{
  int i, j;
  double d1, d2, sx1, cx1, sy1, cy1;
  double snx[7], csx[7], sny[7], csy[7];

  d1 = M_PI*y[0];
  d2 = M_PI*y[1];
  sx1 = sin(d1);
  cx1 = cos(d1);
  sy1 = sin(d2);
  cy1 = cos(d2);
  snx[0] = sx1;
  csx[0] = cx1;
  sny[0] = sy1;
  csy[0] = cy1;
  for (i = 0; i < 6; i++) {
    snx[i + 1] = snx[i] * cx1 + csx[i] * sx1;
    csx[i + 1] = csx[i] * cx1 - snx[i] * sx1;
    sny[i + 1] = sny[i] * cy1 + csy[i] * sy1;
    csy[i + 1] = csy[i] * cy1 - sny[i] * sy1;
  }
  d1 = 0;
  d2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      d1 = d1 + options[index].af[i][j] * snx[i] * sny[j] + options[index].bf[i][j] * csx[i] * csy[j];
      d2 = d2 + options[index].cf[i][j] * snx[i] * sny[j] - options[index].df[i][j] * csx[i] * csy[j];
    }

  // std::ofstream ofstr("data.txt");
  // for (int k = 0; k < 3; k++) {
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "A" << k << "[" << 7 * i + (j + 1) << "] = " << options[k].af[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "B" << k << "[" << 7 * i + (j + 1) << "] = " << options[k].bf[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "C" << k << "[" << 7 * i + (j + 1) << "] = " << options[k].cf[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "D" << k << "[" << 7 * i + (j + 1) << "] = " << options[k].df[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // }
  // ofstr << std::endl;
  // ofstr.close();
  
  return(-sqrt(d1*d1 + d2*d2));
}

// ------------------------------------------------------------------------------------------------
double GrishaginConstrainedProblem::CalculateXDerivative(int index, const vector<double>& y) const
{
  int i, j;
  double dd, d1, d2, t1, t2, sx1, cx1, sy1, cy1;
  double snx[7], csx[7], sny[7], csy[7];

  d1 = M_PI*y[0];
  d2 = M_PI*y[1];
  sx1 = sin(d1);
  cx1 = cos(d1);
  sy1 = sin(d2);
  cy1 = cos(d2);
  snx[0] = sx1;
  csx[0] = cx1;
  sny[0] = sy1;
  csy[0] = cy1;
  for (i = 0; i < 6; i++) {
    snx[i + 1] = snx[i] * cx1 + csx[i] * sx1;
    csx[i + 1] = csx[i] * cx1 - snx[i] * sx1;
    sny[i + 1] = sny[i] * cy1 + csy[i] * sy1;
    csy[i + 1] = csy[i] * cy1 - sny[i] * sy1;
  }
  d1 = 0;
  d2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      d1 = d1 + options[index].af[i][j] * snx[i] * sny[j] + options[index].bf[i][j] * csx[i] * csy[j];
      d2 = d2 + options[index].cf[i][j] * snx[i] * sny[j] - options[index].df[i][j] * csx[i] * csy[j];
    }
  dd = sqrt(d1*d1 + d2*d2);
  t1 = 0;
  t2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      t1 += options[index].af[i][j] * M_PI*(i + 1)*csx[i] * sny[j];
      t1 -= options[index].bf[i][j] * M_PI*(i + 1)*snx[i] * csy[j];
      t2 += options[index].cf[i][j] * M_PI*(i + 1)*csx[i] * sny[j];
      t2 += options[index].df[i][j] * M_PI*(i + 1)*snx[i] * csy[j];
    }
  return(-(t1*d1 + t2*d2) / dd);
}

// ------------------------------------------------------------------------------------------------
double GrishaginConstrainedProblem::CalculateYDerivative(int index, const vector<double>& y) const
{
  int i, j;
  double dd, d1, d2, t1, t2, sx1, cx1, sy1, cy1;
  double snx[7], csx[7], sny[7], csy[7];

  d1 = M_PI*y[0];
  d2 = M_PI*y[1];
  sx1 = sin(d1);
  cx1 = cos(d1);
  sy1 = sin(d2);
  cy1 = cos(d2);
  snx[0] = sx1;
  csx[0] = cx1;
  sny[0] = sy1;
  csy[0] = cy1;
  for (i = 0; i < 6; i++) {
    snx[i + 1] = snx[i] * cx1 + csx[i] * sx1;
    csx[i + 1] = csx[i] * cx1 - snx[i] * sx1;
    sny[i + 1] = sny[i] * cy1 + csy[i] * sy1;
    csy[i + 1] = csy[i] * cy1 - sny[i] * sy1;
  }
  d1 = 0;
  d2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      d1 = d1 + options[index].af[i][j] * snx[i] * sny[j] + options[index].bf[i][j] * csx[i] * csy[j];
      d2 = d2 + options[index].cf[i][j] * snx[i] * sny[j] - options[index].df[i][j] * csx[i] * csy[j];
    }
  dd = sqrt(d1*d1 + d2*d2);
  t1 = 0;
  t2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      t1 += options[index].af[i][j] * M_PI*(j + 1)*snx[i] * csy[j];
      t1 -= options[index].bf[i][j] * M_PI*(j + 1)*csx[i] * sny[j];
      t2 += options[index].cf[i][j] * M_PI*(j + 1)*snx[i] * csy[j];
      t2 += options[index].df[i][j] * M_PI*(j + 1)*csx[i] * sny[j];
    }
  return(-(t1*d1 + t2*d2) / dd);
}

// ------------------------------------------------------------------------------------------------
vector<double> GrishaginConstrainedProblem::ComputeDerivatives(int index, const vector<double>& y) const
{
  vector<double> res(2);
  res[0] = CalculateXDerivative(index, y);
  res[1] = CalculateYDerivative(index, y);
  return res;
}

// ------------------------------------------------------------------------------------------------
double GrishaginConstrainedProblem::rndm20(unsigned char k[])
{
  int i;
  unsigned char k1[45];
  double de2, rndm;

  for (i = 0; i < 38; i++)
    k1[i] = k[i + 7];
  for (i = 38; i < 45; i++)
    k1[i] = 0;
  for (i = 0; i < 45; i++)
    k[i] = (unsigned char)std::abs(k[i] - k1[i]);
  for (i = 27; i < 45; i++)
    k1[i] = k[i - 27];
  for (i = 0; i < 27; i++)
    k1[i] = 0;

  gen(k, k1, 9, 44);
  gen(k, k1, 0, 8);

  rndm = 0.;
  de2 = 1.;
  for (i = 0; i < 36; i++) {
    de2 = de2 / 2;
    rndm = rndm + k[i + 9] * de2;
  }
  return (rndm);
}

// ------------------------------------------------------------------------------------------------
void GrishaginConstrainedProblem::gen(unsigned char k[], unsigned char k1[], int kap1, int kap2)
{
  int jct, i, j;

  jct = 0;
  for (i = kap2; i >= kap1; i--) {
    j = (k[i] + k1[i] + jct) / 2;
    k[i] = k[i] + k1[i] + (unsigned char)jct - (unsigned char)j * 2;
    jct = j;
  }
  if (jct != 0)
    for (i = kap2; i >= kap1; i--) {
      j = (k[i] + jct) / 2;
      k[i] = k[i] + (unsigned char)jct - (unsigned char)j * 2;
      jct = j;
    }
}
