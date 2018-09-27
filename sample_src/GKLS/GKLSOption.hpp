#pragma once
#include <vector>
#include <cmath>
#include <cstdlib>
#include "gkls_random.hpp"

/* Penalty value of the generated function if x is not in D */
#define GKLS_MAX_VALUE        1E+100

/* Value of the machine zero in the floating-point arithmetic */
#define GKLS_PRECISION        1.0E-10

/* Default value of the paraboloid minimum */
#define GKLS_PARABOLOID_MIN   0.0

/* Global minimum value: to be less than GKLS_PARABOLOID_MIN */
#define GKLS_GLOBAL_MIN_VALUE -1.0

/* Max value of the parameter delta for the D2-type class function        */
/* The parameter delta is chosen randomly from (0, GKLS_DELTA_MAX_VALUE ) */
#define GKLS_DELTA_MAX_VALUE  10.0

/* Constant pi */
#ifndef PI
#define PI 3.14159265
#endif

/* Error codes */
#define GKLS_OK                              0
#define GKLS_DIM_ERROR                       1
#define GKLS_NUM_MINIMA_ERROR                2
#define GKLS_FUNC_NUMBER_ERROR               3
#define GKLS_BOUNDARY_ERROR                  4
#define GKLS_GLOBAL_MIN_VALUE_ERROR          5
#define GKLS_GLOBAL_DIST_ERROR               6
#define GKLS_GLOBAL_RADIUS_ERROR             7
#define GKLS_MEMORY_ERROR                    8
#define GKLS_DERIV_EVAL_ERROR                9

/* Reserved error codes */
#define GKLS_GREAT_DIM                      10
#define GKLS_RHO_ERROR                      11
#define GKLS_PEAK_ERROR                     12
#define GKLS_GLOBAL_BASIN_INTERSECTION      13

/* Internal error codes */
#define GKLS_PARABOLA_MIN_COINCIDENCE_ERROR 14
#define GKLS_LOCAL_MIN_COINCIDENCE_ERROR    15
#define GKLS_FLOATING_POINT_ERROR           16

typedef struct {
  double **local_min; /* list of local minimizers coordinates   */
  double *f;          /* list of local minima values            */
  double *w_rho;      /* list of radius weights                 */
  double *peak;       /* list of parameters gamma(i) =          */
                      /*  = local minimum value(i) - paraboloid */
                      /*    minimum within attraction regions   */
                      /*    of local minimizer(i)               */
  double *rho;        /* list of attraction regions radii       */
} T_GKLS_Minima;

/* The structure of type T_GKLS_GlobalMinima contains information      */
/* about the number of global minimizers and their                     */
/* indices in the set of local minimizers                              */

typedef struct {
  unsigned int num_global_minima; /* number of global minima    */
  unsigned int *gm_index;  /* list of indices of generated      */
                           /* minimizers, which are the global ones (elements from 0     */
                           /* to (num_global_minima - 1) of the list) and the local ones */
                           /* (the resting elements of the list)                         */
} T_GKLS_GlobalMinima;

enum GKLSClass { Hard, Simple };
enum GKLSFuncionType { TND, TD, TD2 };

struct GKLSParameters
{
  unsigned dimension;
  double globalMinimumValue;
  unsigned numberOfLocalMinima;
  double globalDistance;
  double globalRadius;
  GKLSFuncionType type;

  GKLSParameters() {}
  GKLSParameters(unsigned _dimension, double _globalMinimumValue,
    unsigned _numberOfLocalMinima, double _globalDistance,
    double _globalRadius, GKLSFuncionType _type) :
    dimension(_dimension),
    globalMinimumValue(_globalMinimumValue),
    numberOfLocalMinima(_numberOfLocalMinima),
    globalDistance(_globalDistance),
    globalRadius(_globalRadius),
    type(_type)
  {}
};

struct GKLSOption
{
  int index;
  bool mIsGeneratorMemoryAllocated;
  bool mIsDomainMemeoryAllocated;
  GKLSFuncionType mFunctionType;
  GKLSRandomGenerator mRndGenerator;

  /*-------------- Variables accessible by the user --------------------- */


  unsigned int GKLS_num_minima; /* number of local minima, >=2  */

  double GKLS_global_dist;  /* distance from the paraboloid minimizer  */
                            /* to the global minimizer                 */
  double GKLS_global_radius;/* radius of the global minimizer          */
                            /* attraction region                       */
  double GKLS_global_value; /* global minimum value,                   */
                            /* test_global_value < GKLS_PARABOLOID_MIN */
  T_GKLS_Minima GKLS_minima;
  /* see the structures type description     */
  T_GKLS_GlobalMinima GKLS_glob;

  /*--------------------------- Global variables ----------------------*/
  int isArgSet; /* isArgSet == 1 if all necessary parameters are set */

  double delta; /* parameter using in D2-type function generation;     */
                /* it is chosen randomly from the                      */
                /* open interval (0,GKLS_DELTA_MAX_VALUE)              */
  unsigned long rnd_counter; /* index of random array elements */

  double* rnd_num, *rand_condition;


  void SetFunctionClass(GKLSClass type, unsigned classDimension)
  {
    GKLS_num_minima = 10;

    if (type == Simple)
    {
      switch (classDimension)
      {
      case 2:
        GKLS_global_dist = 0.9;
        GKLS_global_radius = 0.2;
        break;
      case 3:
        GKLS_global_dist = 0.66;
        GKLS_global_radius = 0.2;
        break;
      case 4:
        GKLS_global_dist = 0.66;
        GKLS_global_radius = 0.2;
        break;
      case 5:
        GKLS_global_dist = 0.66;
        GKLS_global_radius = 0.3;
        break;
      default:
        GKLS_global_dist = 0.66;
        GKLS_global_radius = 0.3;
      }
    }
    else
    {
      switch (classDimension)
      {
      case 2:
        GKLS_global_dist = 0.9;
        GKLS_global_radius = 0.1;
        break;
      case 3:
        GKLS_global_dist = 0.90;
        GKLS_global_radius = 0.2;
        break;
      case 4:
        GKLS_global_dist = 0.90;
        GKLS_global_radius = 0.2;
        break;
      case 5:
        GKLS_global_dist = 0.66;
        GKLS_global_radius = 0.2;
        break;
      default:
        GKLS_global_dist = 0.66;
        GKLS_global_radius = 0.2;
      }
    }
    GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
  }


  void GKLS_free()
  {
    unsigned int i;

    for (i = 0; i < GKLS_num_minima; i++) {
      free(GKLS_minima.local_min[i]);
    }
    free(GKLS_minima.local_min);
    free(GKLS_minima.w_rho);
    free(GKLS_minima.peak);
    free(GKLS_minima.rho);
    free(GKLS_minima.f);
    free(GKLS_glob.gm_index);

    isArgSet = 0;
    mIsGeneratorMemoryAllocated = false;
  }

  int GKLS_alloc(unsigned mDimension)
  {
    unsigned int i;

    if ((mDimension <= 1) || (mDimension >= NUM_RND))
      return GKLS_DIM_ERROR; /* problem dimension error */
    if (GKLS_num_minima <= 1)
      return GKLS_NUM_MINIMA_ERROR; /* erroneous number of local minima */
    if ((GKLS_minima.local_min = (double **)(malloc((size_t)GKLS_num_minima * sizeof(double *)))) == NULL)
      return GKLS_MEMORY_ERROR; /* memory allocation error */
    for (i = 0; i < GKLS_num_minima; i++)
      if ((GKLS_minima.local_min[i] = (double *)(malloc((size_t)mDimension * sizeof(double)))) == NULL)
        return GKLS_MEMORY_ERROR; /* memory allocation error */
    if ((GKLS_minima.w_rho = (double *)(malloc((size_t)GKLS_num_minima * sizeof(double)))) == NULL)
      return GKLS_MEMORY_ERROR;   /* memory allocation error */
    if ((GKLS_minima.peak = (double *)(malloc((size_t)GKLS_num_minima * sizeof(double)))) == NULL)
      return GKLS_MEMORY_ERROR;   /* memory allocation error */
    if ((GKLS_minima.rho = (double *)(malloc((size_t)GKLS_num_minima * sizeof(double)))) == NULL)
      return GKLS_MEMORY_ERROR;   /* memory allocation error */
    if ((GKLS_minima.f = (double *)(malloc((size_t)GKLS_num_minima * sizeof(double)))) == NULL)
      return GKLS_MEMORY_ERROR;   /* memory allocation error */
    if ((GKLS_glob.gm_index =
      (unsigned int *)(malloc((size_t)GKLS_num_minima * sizeof(unsigned int)))) == NULL)
      return GKLS_MEMORY_ERROR;   /* memory allocation error */
    else
      GKLS_glob.num_global_minima = 0;
    mIsGeneratorMemoryAllocated = true;
    return GKLS_OK; /* no errors */
  }


  int GKLS_parameters_check(unsigned mDimension, std::vector<double> mLoBound, std::vector<double> mUpBound) const
  {
    unsigned int i;
    double min_side, tmp;

    if (GKLS_num_minima <= 1)  /* number of local minima error */
      return GKLS_NUM_MINIMA_ERROR;


    if (GKLS_global_value >= GKLS_PARABOLOID_MIN - GKLS_PRECISION)
      return GKLS_GLOBAL_MIN_VALUE_ERROR; /* the global minimum value must   */
                                          /* be less than the paraboloid min */
                                          /* Find min_side = min |b(i)-a(i)|, D=[a,b], and                   */
                                          /* check the distance from paraboloid vertex to global minimizer   */
    min_side = mUpBound[0] - mLoBound[0];
    for (i = 1; i < mDimension; i++)
      if ((tmp = mUpBound[i] - mLoBound[i]) < min_side)
        min_side = tmp;
    if ((GKLS_global_dist >= 0.5*min_side - GKLS_PRECISION) ||
      (GKLS_global_dist <= GKLS_PRECISION))
      return GKLS_GLOBAL_DIST_ERROR; /* global distance error */
    if ((GKLS_global_radius >= 0.5*GKLS_global_dist + GKLS_PRECISION) ||
      (GKLS_global_radius <= GKLS_PRECISION))
      return GKLS_GLOBAL_RADIUS_ERROR; /* global minimizer attr. radius error */

    return GKLS_OK; /* no errors */
  }

  int GKLS_arg_generate(unsigned int nf, double* argmin, unsigned mDimension, std::vector<double> mLoBound, std::vector<double> mUpBound)
  {
    unsigned int i, j;
    int error;
    double sin_phi; /* for generating of the global minimizer coordinates */
                    /* by using the generalized spherical coordinates     */
    double gap = GKLS_global_radius; /* gap > 0 */
                                            /* the minimal distance of any local minimizer to the attraction*/
                                            /* region of the global minimizer M(1); the value               */
                                            /* GKLS_global_radius is given by default and can be changed,   */
                                            /* but it should not be too small.                              */

                                            /* Check function number */
    if ((nf < 1) || (nf > 100)) return GKLS_FUNC_NUMBER_ERROR;

    /* Check parameters */
    if ((error = GKLS_parameters_check(mDimension, mLoBound, mUpBound)) != GKLS_OK) return error;

    /* Allocate memory */
    if ((error = GKLS_alloc(mDimension)) != GKLS_OK) return error;

    /* Set random seed */
    if ((error = GKLS_initialize_rnd(mDimension, GKLS_num_minima, nf)) != GKLS_OK)
      return error;

    mRndGenerator.GenerateNextNumbers();//ranf_array(rnd_num, NUM_RND); /* get random sequence */
    rnd_counter = 0L;   /* index of the random element from the sequence */
                               /* to be used as the next random number          */

                               /* Set the paraboloid minimizer coordinates and */
                               /* the paraboloid minimum value                 */
    for (i = 0; i < mDimension; i++) {
      GKLS_minima.local_min[0][i] = mLoBound[i] +
        rnd_num[rnd_counter] * (mUpBound[i] - mLoBound[i]);
      rnd_counter++;
      if (rnd_counter == NUM_RND) {
        mRndGenerator.GenerateNextNumbers();//ranf_array(rnd_num, NUM_RND);
        rnd_counter = 0L;
      }
    } /* for coordinates */
    GKLS_minima.f[0] = GKLS_PARABOLOID_MIN; /* fix the paraboloid min value */

                                                   /* Generate the global minimizer using generalized spherical coordinates*/
                                                   /* with an arbitrary vector phi and the fixed radius GKLS_global_radius */

                                                   /* First, generate an angle 0 <= phi(0) <= PI, and the coordinate x(0)*/
    mRndGenerator.GenerateNextNumbers(); //ranf_array(rnd_num, NUM_RND);
    rnd_counter = 0L;

    if (argmin != 0)
    {
      for (j = 0; j < mDimension; j++)
        GKLS_minima.local_min[1][j] = argmin[j];

      sin_phi = sin(PI*rnd_num[rnd_counter]);
      rnd_counter++;

      for (j = 1; j < mDimension - 1; j++) {
        sin_phi *= sin(2.0*PI*rnd_num[rnd_counter]);
        rnd_counter++;
      }
    }
    else
    {
      GKLS_minima.local_min[1][0] = GKLS_minima.local_min[0][0] +
        GKLS_global_dist*cos(PI*rnd_num[rnd_counter]);
      if ((GKLS_minima.local_min[1][0] > mUpBound[0] - GKLS_PRECISION) ||
        (GKLS_minima.local_min[1][0] < mLoBound[0] + GKLS_PRECISION))
        GKLS_minima.local_min[1][0] = GKLS_minima.local_min[0][0] -
        GKLS_global_dist*cos(PI*rnd_num[rnd_counter]);
      sin_phi = sin(PI*rnd_num[rnd_counter]);
      rnd_counter++;

      /* Generate the remaining angles 0<=phi(i)<=2*PI, and         */
      /* the coordinates x(i), i=1,...,mDimension-2 (not last!)       */
      for (j = 1; j < mDimension - 1; j++) {
        GKLS_minima.local_min[1][j] = GKLS_minima.local_min[0][j] +
          GKLS_global_dist*cos(2.0*PI*rnd_num[rnd_counter])*sin_phi;
        if ((GKLS_minima.local_min[1][j] > mUpBound[j] - GKLS_PRECISION) ||
          (GKLS_minima.local_min[1][j] < mLoBound[j] + GKLS_PRECISION))
          GKLS_minima.local_min[1][j] = GKLS_minima.local_min[0][j] -
          GKLS_global_dist*cos(2.0*PI*rnd_num[rnd_counter])*sin_phi;
        sin_phi *= sin(2.0*PI*rnd_num[rnd_counter]);
        rnd_counter++;
      }

      /* Generate the last coordinate x(mDimension-1) */
      GKLS_minima.local_min[1][mDimension - 1] = GKLS_minima.local_min[0][mDimension - 1] +
        GKLS_global_dist*sin_phi;
      if ((GKLS_minima.local_min[1][mDimension - 1] > mUpBound[mDimension - 1] - GKLS_PRECISION) ||
        (GKLS_minima.local_min[1][mDimension - 1] < mLoBound[mDimension - 1] + GKLS_PRECISION))
        GKLS_minima.local_min[1][mDimension - 1] =
        GKLS_minima.local_min[0][mDimension - 1] - GKLS_global_dist*sin_phi;
    }
    /* Set the global minimum value */
    GKLS_minima.f[1] = GKLS_global_value;

    /* Set the weight coefficients w_rho(i) */
    for (i = 0; i < GKLS_num_minima; i++)
      GKLS_minima.w_rho[i] = 0.99;
    GKLS_minima.w_rho[1] = 1.0;


    /* Set the parameter delta for D2-type functions       */
    /* It is chosen randomly from (0,GKLS_DELTA_MAX_VALUE) */
    delta = GKLS_DELTA_MAX_VALUE*rnd_num[rnd_counter];
    rnd_counter++;
    if (rnd_counter == NUM_RND) {
      mRndGenerator.GenerateNextNumbers(); //ranf_array(rnd_num, NUM_RND);
      rnd_counter = 0L;
    }

    /* Choose randomly coordinates of local minimizers       */
    /* This procedure is repeated while the local minimizers */
    /* coincide (external do...while);                       */
    /* The internal cycle do..while serves to choose local   */
    /* minimizers in certain distance from the attraction    */
    /* region of the global minimizer M(i)                   */
    do
    {
      i = 2;
      while (i < GKLS_num_minima) {
        do
        {
          mRndGenerator.GenerateNextNumbers();//ranf_array(rnd_num, NUM_RND);
          rnd_counter = 0L;
          for (j = 0; j < mDimension; j++) {
            GKLS_minima.local_min[i][j] = mLoBound[j] +
              rnd_num[rnd_counter] * (mUpBound[j] - mLoBound[j]);
            rnd_counter++;
            if (rnd_counter == NUM_RND) {
              mRndGenerator.GenerateNextNumbers(); //ranf_array(rnd_num, NUM_RND);
              rnd_counter = 0L;
            }
          }
          /* Check wether this local minimizer belongs to a zone of */
          /* the global minimizer M(i)                              */
        } while ((GKLS_global_radius + gap) -
          GKLS_norm(GKLS_minima.local_min[i], GKLS_minima.local_min[1], mDimension)
          > GKLS_PRECISION);
        i++;
      }
      error = GKLS_coincidence_check(mDimension);
    } while ((error == GKLS_PARABOLA_MIN_COINCIDENCE_ERROR) ||
      (error == GKLS_LOCAL_MIN_COINCIDENCE_ERROR));
    error = GKLS_set_basins(mDimension, mLoBound, mUpBound);
    if (error == GKLS_OK) isArgSet = 1; /* All the parameters are set */
                                               /* and the user can Calculate a specific test function or          */
                                               /* its partial derivative by calling corresponding subroutines    */

    return error;
  }


  int GKLS_coincidence_check(unsigned mDimension) const
  {
    unsigned int i, j;

    /* Check wether some local minimizer coincides with the paraboloid minimizer */
    for (i = 2; i < GKLS_num_minima; i++)
    {
      if (GKLS_norm(GKLS_minima.local_min[i], GKLS_minima.local_min[0], mDimension) < GKLS_PRECISION)
        return GKLS_PARABOLA_MIN_COINCIDENCE_ERROR;
    }

    /* Check wether there is a pair of identical local minimizers */
    for (i = 1; i < GKLS_num_minima - 1; i++)
      for (j = i + 1; j < GKLS_num_minima; j++)
      {
        if (GKLS_norm(GKLS_minima.local_min[i], GKLS_minima.local_min[j], mDimension) < GKLS_PRECISION)
          return GKLS_LOCAL_MIN_COINCIDENCE_ERROR;
      }

    return GKLS_OK;
  }

  int GKLS_set_basins(unsigned mDimension, std::vector<double> mLoBound, std::vector<double> mUpBound)
  {
    unsigned int i, j;
    double temp_min;         /*  temporary  */
    double temp_d1, temp_d2; /*  variables  */
    double dist;  /* for finding the distance between two minimizers */

                  /****************************************************************************/
                  /* First, set the radii rho(i) of the attraction regions: these values are  */
                  /* defined in such a way that attraction regions are as large as possible   */
                  /* and do not overlap; it is not required that the attraction region of each*/
                  /* local minimizer be entirely contained in D. The values found in such     */
                  /* a manner are corrected then by the weight coefficients w(i)              */
                  /****************************************************************************/

                  /* Calculate dist(i) - the minimal distance from the minimizer i to         */
                  /*                     the other minimizers.                                */
                  /* Set the initial value of rho(i) as rho(i) = dist(i)/2: so the attraction */
                  /* regions do not overlap                                                   */
    for (i = 0; i < GKLS_num_minima; i++)
    {
      temp_min = GKLS_MAX_VALUE;
      for (j = 0; j < GKLS_num_minima; j++)
        if (i != j)
        {
          if ((temp_d1 = GKLS_norm(GKLS_minima.local_min[i], GKLS_minima.local_min[j], mDimension)) < temp_min)
            temp_min = temp_d1;
        }

      dist = temp_min / 2.0;

      GKLS_minima.rho[i] = dist;
    }

    /* Since the radius of the attraction region of the global minimizer            */
    /* is fixed by the user, the generator adjusts the radii of the attraction      */
    /* regions, eventually overlapping with the attraction region of the global     */
    /* minimizer. To do this, it checks whether the attraction region radius        */
    /* of each local minimizer exceeds the distance between this minimizer          */
    /* and the attraction region of the global minimizer.                           */
    /* If such a situation is verified the generator decreases the attraction       */
    /* region radius of the local minimizer setting it equal to he distance between */
    /* the local minimizer and the attraction region of the global minimizer.       */
    /* Note that the radius of the attraction region of the global minimizer can    */
    /* not be greater than one half of the distance (defined by the user) between   */
    /* the global minimizer and the paraboloid vertex. So, the initially defined    */
    /* attraction regions of the global minimizer and the paraboloid vertex do not  */
    /* overlap even when the global minimizer is the closest minimizer to the       */
    /* paraboloid vertex.                                                           */
    GKLS_minima.rho[1] = GKLS_global_radius;
    for (i = 2; i < GKLS_num_minima; i++)
    {
      if ((dist = (GKLS_norm(GKLS_minima.local_min[i], GKLS_minima.local_min[1], mDimension)
        - GKLS_global_radius - GKLS_PRECISION)) < GKLS_minima.rho[i])
        GKLS_minima.rho[i] = dist;
    }
    /* Try to expand the attraction regions of local minimizers until they      */
    /* do not overlap                                                           */
    for (i = 0; i < GKLS_num_minima; i++)
    {
      if (i != 1) /* The radius of the attr. region of the global min is fixed  */
      { /* rho(i) := max {rho(i),min[||M(i)-M(j)|| - rho(j): i !=j]},      */
        temp_min = GKLS_MAX_VALUE;
        for (j = 0; j < GKLS_num_minima; j++)
          if (i != j)
          {
            if ((temp_d1 = GKLS_norm(GKLS_minima.local_min[i], GKLS_minima.local_min[j], mDimension) - GKLS_minima.rho[j]) < temp_min)
              temp_min = temp_d1;
          }
        /* Increase the radius rho(i) if it is possible */
        if (temp_min > GKLS_minima.rho[i] + GKLS_PRECISION)
          GKLS_minima.rho[i] = temp_min;
      }
    }

    /* Correct the radii by weight coefficients w(i)                    */
    /* The weight coefficients can be chosen randomly;                  */
    /* here they are defined by default as:                             */
    /*    w(i) = 0.99, i != 1 , and w(1) = 1.0 (global min index = 1)   */
    for (i = 0; i < GKLS_num_minima; i++)
      GKLS_minima.rho[i] = GKLS_minima.w_rho[i] * GKLS_minima.rho[i];

    /*********************************************************************/
    /* Set the local minima values f(i) of test functions as follows:    */
    /*   f(i) = cond_min(i) - peak(i), i != 1 (global min index = 1)     */
    /*   f(0) = GKLS_PARABOLOID_MIN, f(1) = GKLS_GLOBAL_MIN_VALUE,       */
    /* where cond_min(i) is the paraboloid minimum value at the boundary */
    /* B={||x-M(i)||=rho(i)} of the attraction region of the local       */
    /* minimizer M(i), i.e.                                              */
    /*  cond_min(i) =                                                    */
    /*  = paraboloid g() value at (M(i)+rho(i)*(T-M(i))/norm(T-M)) =     */
    /*  = (rho(i) - norm(T-M(i)))^2 + t,                                 */
    /*  g(x) = ||x-T||^2 + t, x in D of R^mDimension                       */
    /*                                                                   */
    /*  The values of peak(i) are chosen randomly from an interval       */
    /* (0, 2*rho(i), so that the values f(i) depend on radii rho(i) of   */
    /* the attraction regions, 2<=i<mDimension.                            */
    /*  The condition f(x*)=f(1) <= f(i) must be satisfied               */
    /*********************************************************************/
    /* Fix to 0 the peak(i) values of the paraboloid and the global min  */
    GKLS_minima.peak[0] = 0.0; GKLS_minima.peak[1] = 0.0;
    for (i = 2; i < GKLS_num_minima; i++) {
      /* Set values peak(i), i>= 2 */
      /* Note that peak(i) is such that the function value f(i) is smaller*/
      /* than min(GKLS_GLOBAL_MIN_VALUE, 2*rho(i))                        */
      temp_d1 = GKLS_norm(GKLS_minima.local_min[0], GKLS_minima.local_min[i], mDimension);
      temp_min = (GKLS_minima.rho[i] - temp_d1)*(GKLS_minima.rho[i] - temp_d1) +
        GKLS_minima.f[0]; /*the conditional minimum at the boundary*/

      temp_d1 = (1.0 + rnd_num[rnd_counter])*GKLS_minima.rho[i];
      temp_d2 = rnd_num[rnd_counter] * (temp_min - GKLS_global_value);
      /* temp_d1 := min(temp_d1, temp_d2) */
      if (temp_d2 < temp_d1) temp_d1 = temp_d2;
      GKLS_minima.peak[i] = temp_d1;

      rnd_counter++;
      if (rnd_counter == NUM_RND)
      {
        mRndGenerator.GenerateNextNumbers(); //ranf_array(rnd_num, NUM_RND);
        rnd_counter = 0L;
      }

      GKLS_minima.f[i] = temp_min - GKLS_minima.peak[i];
    }

    /*********************************************************************/
    /* Find all possible global minimizers and                           */
    /* create a list of their indices among all the minimizers           */
    /* Note that the paraboloid minimum can not be the global one because*/
    /* the global optimum value is set to be less than the paraboloid    */
    /* minimum value                                                     */
    /*********************************************************************/
    GKLS_glob.num_global_minima = 0;
    for (i = 0; i < GKLS_num_minima; i++)
      if ((GKLS_minima.f[i] >= GKLS_global_value - GKLS_PRECISION) &&
        (GKLS_minima.f[i] <= GKLS_global_value + GKLS_PRECISION))
      {
        GKLS_glob.gm_index[GKLS_glob.num_global_minima] = i;
        GKLS_glob.num_global_minima++;
        /* The first GKLS_glob.num_global_minima elements of the list    */
        /* contain the indices of the global minimizers                  */
      }
      else
        GKLS_glob.gm_index[GKLS_num_minima - 1 - i + GKLS_glob.num_global_minima]
        = i;
    /* The remaining elements of the list                            */
    /* contain the indices of local (non global) minimizers          */

    if (GKLS_glob.num_global_minima == 0) /*   erroneous case:       */
      return GKLS_FLOATING_POINT_ERROR;  /* some programmer's error */

    return GKLS_OK;
  }
  int GKLS_initialize_rnd(unsigned int dim, unsigned int nmin, int nf)
  {
    long seed;
    /* seed number between 0 and 2^30-3 = 1,073,741,821*/

    seed = (nf - 1) + (nmin - 1) * 100 + dim * 1000000L;
    /* If big values of nmin and dim are required, */
    /* one must check wether seed <= 1073741821    */

    mRndGenerator.Initialize(seed, rnd_num, rand_condition);//ranf_start(seed);

    return GKLS_OK;
  }
  double GKLS_norm(const double *x1, const double *x2, unsigned mDimension) const
  {
    unsigned int i;
    double norm = 0.0;
    for (i = 0; i < mDimension; i++)
      norm += (x1[i] - x2[i])*(x1[i] - x2[i]);
    return sqrt(norm);
  }
};
