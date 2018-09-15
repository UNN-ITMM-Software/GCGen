#pragma once

/// Grishagin problems parameters
struct GrishaginOption
{
  unsigned char icnf[45];
  double af[7][7], bf[7][7], cf[7][7], df[7][7];
  int index;
};

extern const double rand_minimums[];
extern const unsigned char matcon[10][45];
