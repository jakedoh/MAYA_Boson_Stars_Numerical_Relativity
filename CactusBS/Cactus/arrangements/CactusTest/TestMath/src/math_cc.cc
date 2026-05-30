#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Arguments.h>



// Test long double only if it is actually used, i.e. if e.g.
// CCTK_REAL16 is long double
#if defined SIZEOF_LONG_DOUBLE && HAVE_CCTK_REAL16 && SIZEOF_LONG_DOUBLE==16
#  define TEST_LONG_DOUBLE
#endif



extern "C"
void TestMath_CC(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
#ifdef TEST_LONG_DOUBLE
  // Ensure that long double works in C++
  int check_long_double[sizeof(long double) == SIZEOF_LONG_DOUBLE ? +1 : -1]
    CCTK_ATTRIBUTE_UNUSED;
#endif

  int i CCTK_ATTRIBUTE_UNUSED = 1;
  float f CCTK_ATTRIBUTE_UNUSED = 2;
  double d CCTK_ATTRIBUTE_UNUSED = 3;
#ifdef TEST_LONG_DOUBLE
  long double ld CCTK_ATTRIBUTE_UNUSED = 4;
#endif
  
  // Note: In C, copysign and signbit are functions with different
  // names for each type. In C++, they are overloaded, and have the
  // same name for each type.
  
  // Test calling the functions with a "std::" prefix
  
#ifdef HAVE_COPYSIGN
  f = std::copysign(1.0f, 1.0f);
  d = std::copysign(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::copysign(1.0L, 1.0L);
#  endif
#endif
  
#if 0                           // does not exist
  f = std::fabs(1.0f);
  d = std::fabs(1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::fabs(1.0L);
#  endif
#endif
  
  // Ensure that abs does not return an integer type
  int check_std_absf[std::abs(0.1f) > 0 ? +1 : -1] CCTK_ATTRIBUTE_UNUSED;
  f = std::abs(1.0f);
  int check_std_abs[std::abs(0.1) > 0 ? +1 : -1] CCTK_ATTRIBUTE_UNUSED;
  d = std::abs(1.0);
#  ifdef TEST_LONG_DOUBLE
  int check_std_absl[std::abs(0.1L) > 0 ? +1 : -1] CCTK_ATTRIBUTE_UNUSED;
  ld = std::abs(1.0L);
#  endif
  
#if 0                           // does not exist
  f = std::fmax(1.0f, 1.0f);
  d = std::fmax(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::fmax(1.0L, 1.0L);
#  endif
#endif
  
  f = std::max(1.0f, 1.0f);
  d = std::max(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::max(1.0L, 1.0L);
#  endif
  
#if 0                           // does not exist
  f = std::fmin(1.0f, 1.0f);
  d = std::fmin(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::fmin(1.0L, 1.0L);
#  endif
#endif
  
  f = std::min(1.0f, 1.0f);
  d = std::min(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::min(1.0L, 1.0L);
#  endif

#ifdef HAVE_FPCLASSIFY
  i = std::fpclassify(1.0f);
  i = std::fpclassify(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = std::fpclassify(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISFINITE
  i = std::isfinite(1.0f);
  i = std::isfinite(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = std::isfinite(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISINF
  i = std::isinf(1.0f);
  i = std::isinf(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = std::isinf(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISNAN
  i = std::isnan(1.0f);
  i = std::isnan(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = std::isnan(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISNORMAL
  i = std::isnormal(1.0f);
  i = std::isnormal(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = std::isnormal(1.0L);
#  endif
#endif
  
  f = std::pow(1.0f, 1.0f);
  d = std::pow(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = std::pow(1.0L, 1.0L);
#  endif
  
#ifdef HAVE_SIGNBIT
  i = std::signbit(1.0f);
  i = std::signbit(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = std::signbit(1.0L);
#  endif
#endif
  
  // Also test "using namespace std", without any prefix
  
  using namespace std;
  
#ifdef HAVE_COPYSIGN
  f = copysign(1.0f, 1.0f);
  d = copysign(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = copysign(1.0L, 1.0L);
#  endif
#endif
  
  f = fabs(1.0f);
  d = fabs(1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = fabs(1.0L);
#  endif
  
  // Ensure that abs does not return an integer type
  int check_absf[abs(0.1f) > 0 ? +1 : -1] CCTK_ATTRIBUTE_UNUSED;
  f = abs(1.0f);
  int check_abs[abs(0.1) > 0 ? +1 : -1] CCTK_ATTRIBUTE_UNUSED;
  d = abs(1.0);
#  ifdef TEST_LONG_DOUBLE
  int check_absl[abs(0.1L) > 0 ? +1 : -1] CCTK_ATTRIBUTE_UNUSED;
  ld = abs(1.0L);
#  endif
  
  f = fmax(1.0f, 1.0f);
  d = fmax(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = fmax(1.0L, 1.0L);
#  endif
  
  f = max(1.0f, 1.0f);
  d = max(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = max(1.0L, 1.0L);
#  endif
  
  f = fmin(1.0f, 1.0f);
  d = fmin(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = fmin(1.0L, 1.0L);
#  endif
  
  f = min(1.0f, 1.0f);
  d = min(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = min(1.0L, 1.0L);
#  endif
  
#ifdef HAVE_FPCLASSIFY
  i = fpclassify(1.0f);
  i = fpclassify(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = fpclassify(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISFINITE
  i = isfinite(1.0f);
  i = isfinite(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = isfinite(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISINF
  i = isinf(1.0f);
  i = isinf(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = isinf(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISNAN
  i = isnan(1.0f);
  i = isnan(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = isnan(1.0L);
#  endif
#endif
  
#ifdef HAVE_ISNORMAL
  i = isnormal(1.0f);
  i = isnormal(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = isnormal(1.0L);
#  endif
#endif
  
  f = pow(1.0f, 1.0f);
  d = pow(1.0, 1.0);
#  ifdef TEST_LONG_DOUBLE
  ld = pow(1.0L, 1.0L);
#  endif
  
#ifdef HAVE_SIGNBIT
  i = signbit(1.0f);
  i = signbit(1.0);
#  ifdef TEST_LONG_DOUBLE
  i = signbit(1.0L);
#  endif
#endif
  
}
