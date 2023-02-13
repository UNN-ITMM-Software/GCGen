#define _USE_MATH_DEFINES

#include <math.h>
#include <algorithm>
#include <cassert>

#include "grishagin_function.hpp"

// #include <fstream>

// ------------------------------------------------------------------------------------------------
TGrishaginProblem::TGrishaginProblem(int problemIndex) : IOptProblem()
{
  mProblemIndex = problemIndex;
  option.index = problemIndex;
  SetOptions();
  mDimension = 2;
  mLoBound = { 0.0, 0.0 };
  mUpBound = { 1.0, 1.0 };
  mOptimumPoint = { rand_minimums[2 * (mProblemIndex - 1)], rand_minimums[2 * (mProblemIndex - 1) + 1] };
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
void TGrishaginProblem::SetOptions()
{
  assert(mProblemIndex > 0 && mProblemIndex <= 100);

  int lst, i, j, i1, i2, i3;
  int nf = mProblemIndex;

  if (nf < 1 || nf>100)
    nf = 1;
  lst = 10;
  i1 = (nf - 1) / lst;
  i2 = i1*lst;
  for (j = 0; j < 45; j++)
    option.icnf[j] = matcon[i1][j];
  if (i2 != (nf - 1)) {
    i3 = nf - 1 - i2;
    for (j = 1; j <= i3; j++)
      for (i = 0; i < 196; i++)
        rndm20(option.icnf);
  }
  for (j = 0; j < 7; j++)
    for (i = 0; i < 7; i++) {
      option.af[i][j] = 2.*rndm20(option.icnf) - 1.;
      option.cf[i][j] = 2.*rndm20(option.icnf) - 1.;
    }
  for (j = 0; j < 7; j++)
    for (i = 0; i < 7; i++) {
      option.bf[i][j] = 2.*rndm20(option.icnf) - 1.;
      option.df[i][j] = 2.*rndm20(option.icnf) - 1.;
    }
}

// ------------------------------------------------------------------------------------------------
double TGrishaginProblem::Compute(int index, const vector<double>& y) const
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
      d1 = d1 + option.af[i][j] * snx[i] * sny[j] + option.bf[i][j] * csx[i] * csy[j];
      d2 = d2 + option.cf[i][j] * snx[i] * sny[j] - option.df[i][j] * csx[i] * csy[j];
    }

  // std::ofstream ofstr("data.txt");
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "A" << "[" << 7 * i + (j + 1) << "] = " << options.af[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "B" << "[" << 7 * i + (j + 1) << "] = " << options.bf[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "C" << "[" << 7 * i + (j + 1) << "] = " << options.cf[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     ofstr << "D" << "[" << 7 * i + (j + 1) << "] = " << options.df[i][j] << std::endl;
  //   }
  // }
  // ofstr << std::endl;
  // ofstr.close();
  
  return(-sqrt(d1*d1 + d2*d2));
}

// ------------------------------------------------------------------------------------------------
double TGrishaginProblem::CalculateXDerivative(const vector<double>& y) const
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
      d1 = d1 + option.af[i][j] * snx[i] * sny[j] + option.bf[i][j] * csx[i] * csy[j];
      d2 = d2 + option.cf[i][j] * snx[i] * sny[j] - option.df[i][j] * csx[i] * csy[j];
    }
  dd = sqrt(d1*d1 + d2*d2);
  t1 = 0;
  t2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      t1 += option.af[i][j] * M_PI*(i + 1)*csx[i] * sny[j];
      t1 -= option.bf[i][j] * M_PI*(i + 1)*snx[i] * csy[j];
      t2 += option.cf[i][j] * M_PI*(i + 1)*csx[i] * sny[j];
      t2 += option.df[i][j] * M_PI*(i + 1)*snx[i] * csy[j];
    }
  return(-(t1*d1 + t2*d2) / dd);
}

// ------------------------------------------------------------------------------------------------
double TGrishaginProblem::CalculateYDerivative(const vector<double>& y) const
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
      d1 = d1 + option.af[i][j] * snx[i] * sny[j] + option.bf[i][j] * csx[i] * csy[j];
      d2 = d2 + option.cf[i][j] * snx[i] * sny[j] - option.df[i][j] * csx[i] * csy[j];
    }
  dd = sqrt(d1*d1 + d2*d2);
  t1 = 0;
  t2 = 0;
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++) {
      t1 += option.af[i][j] * M_PI*(j + 1)*snx[i] * csy[j];
      t1 -= option.bf[i][j] * M_PI*(j + 1)*csx[i] * sny[j];
      t2 += option.cf[i][j] * M_PI*(j + 1)*snx[i] * csy[j];
      t2 += option.df[i][j] * M_PI*(j + 1)*csx[i] * sny[j];
    }
  return(-(t1*d1 + t2*d2) / dd);
}

// ------------------------------------------------------------------------------------------------
vector<double> TGrishaginProblem::ComputeDerivatives(int index, const vector<double>& y) const
{
  vector<double> res(2);
  res[0] = CalculateXDerivative(y);
  res[1] = CalculateYDerivative(y);
  return res;
}

// ------------------------------------------------------------------------------------------------
double TGrishaginProblem::rndm20(unsigned char k[])
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
void TGrishaginProblem::gen(unsigned char k[], unsigned char k1[], int kap1, int kap2)
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
