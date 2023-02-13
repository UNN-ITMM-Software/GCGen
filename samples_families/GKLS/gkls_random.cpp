#include "gkls_random.hpp"



GKLSRandomGenerator::GKLSRandomGenerator() {}

GKLSRandomGenerator::~GKLSRandomGenerator() {}

void GKLSRandomGenerator::Initialize(long seed, double* rnd_num_mem, double* rand_condition_mem)
{
  ran_u = rand_condition_mem;
  rnd_num = rnd_num_mem;
  int j;
  long t, s;
  double u[KK + KK - 1], ul[KK + KK - 1];
  double ulp = (1.0 / (1L << 30)) / (1L << 22);          /* 2 to the -52 */
  double ss = 2.0*ulp*((seed & 0x3fffffff) + 2);

  for (j = 0; j < KK; j++) {
    u[j] = ss; ul[j] = 0.0;                  /* bootstrap the buffer */
    ss += ss; if (ss >= 1.0) ss -= 1.0 - 2 * ulp;  /* cyclic shift of 51 bits */
  }
  for (; j < KK + KK - 1; j++) u[j] = ul[j] = 0.0;
  u[1] += ulp; ul[1] = ulp;         /* make u[1] (and only u[1]) "odd" */
  s = seed & 0x3fffffff;
  t = TT - 1; while (t) {
    for (j = KK - 1; j > 0; j--) ul[j + j] = ul[j], u[j + j] = u[j];   /* "square" */
    for (j = KK + KK - 2; j > KK - LL; j -= 2)
      ul[KK + KK - 1 - j] = 0.0, u[KK + KK - 1 - j] = u[j] - ul[j];
    for (j = KK + KK - 2; j >= KK; j--) if (ul[j]) {
      ul[j - (KK - LL)] = ulp - ul[j - (KK - LL)],
        u[j - (KK - LL)] = mod_sum(u[j - (KK - LL)], u[j]);
      ul[j - KK] = ulp - ul[j - KK], u[j - KK] = mod_sum(u[j - KK], u[j]);
    }
    if (is_odd(s)) {                          /* "multiply by z" */
      for (j = KK; j > 0; j--)  ul[j] = ul[j - 1], u[j] = u[j - 1];
      ul[0] = ul[KK], u[0] = u[KK];       /* shift the buffer cyclically */
      if (ul[KK]) ul[LL] = ulp - ul[LL], u[LL] = mod_sum(u[LL], u[KK]);
    }
    if (s) s >>= 1; else t--;
  }
  for (j = 0; j < LL; j++) ran_u[j + KK - LL] = u[j];
  for (; j < KK; j++) ran_u[j - LL] = u[j];
}

void GKLSRandomGenerator::GenerateNextNumbers()
{
  int i, j, n = NUM_RND;
  for (j = 0; j < KK; j++) rnd_num[j] = ran_u[j];
  for (; j < n; j++) rnd_num[j] = mod_sum(rnd_num[j - KK], rnd_num[j - LL]);

  for (i = 0; i < LL; i++, j++) ran_u[i] = mod_sum(rnd_num[j - KK], rnd_num[j - LL]);
  for (; i < KK; i++, j++) ran_u[i] = mod_sum(rnd_num[j - KK], ran_u[i - LL]);
}

double GKLSRandomGenerator::GetRandomNumber(int indx) const
{
  return rnd_num[indx];
}
double GKLSRandomGenerator::GetGeneratorState(int indx) const
{
  return ran_u[indx];
}
