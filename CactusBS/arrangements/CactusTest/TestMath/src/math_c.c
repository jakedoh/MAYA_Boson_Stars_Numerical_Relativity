#include <math.h>

#include <cctk.h>
#include <cctk_Arguments.h>



void TestMath_C(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  volatile int         i  CCTK_ATTRIBUTE_UNUSED = 1;
  volatile float       f  CCTK_ATTRIBUTE_UNUSED = 2;
  volatile double      d  CCTK_ATTRIBUTE_UNUSED = 3;
#ifdef SIZEOF_LONG_DOUBLE
  volatile long double ld CCTK_ATTRIBUTE_UNUSED = 4;
#endif
  
  /* Note: copysign is a function, with different names for each type.
     isnan and friends are macros, using the same name for each
     type. */
  
#ifdef HAVE_COPYSIGN
  f = copysignf(1.0f, 1.0f);
  d = copysign(1.0, 1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  ld = copysignl(1.0L, 1.0L);
#  endif
#endif
  
  f = fabsf(1.0f);
  d = fabs(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  ld = fabsl(1.0L);
#  endif
  
  f = fmaxf(1.0f, 1.0f);
  d = fmax(1.0, 1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  ld = fmaxl(1.0L, 1.0L);
#  endif
  
  f = fminf(1.0f, 1.0f);
  d = fmin(1.0, 1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  ld = fminl(1.0L, 1.0L);
#  endif
  
#ifdef HAVE_FPCLASSIFY
  i = fpclassify(1.0f);
  i = fpclassify(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  i = fpclassify(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISFINITE
  i = isfinite(1.0f);
  i = isfinite(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  i = isfinite(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISINF
  i = isinf(1.0f);
  i = isinf(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  i = isinf(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISNAN
  i = isnan(1.0f);
  i = isnan(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  i = isnan(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISNORMAL
  i = isnormal(1.0f);
  i = isnormal(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  i = isnormal(1.0L);
#  endif
#endif
  
  f = powf(1.0f, 1.0f);
  d = pow(1.0, 1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  ld = powl(1.0L, 1.0L);
#  endif
  
#ifdef HAVE_SIGNBIT
  i = signbit(1.0f);
  i = signbit(1.0);
#  ifdef SIZEOF_LONG_DOUBLE
  i = signbit(1.0L);
#  endif
#endif
  
}
