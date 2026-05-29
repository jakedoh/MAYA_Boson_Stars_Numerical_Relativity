#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "loopcontrol.h"
#include "vectors.h"
#include "FLRWBackground.h"

#define KRANC_C

namespace FLRWBackground {

/*
Solves the equation eps^3(1+eps) = 27 K^3 u a^-4, where eps is the variable,
to get the intrinsic energy for the hydro equation of state.
Uses standard Newton-Raphson algorithm.
*/

extern "C" double getEps(double u , double K , double a, double tol) {

  double xk, xkplusone, f, fprime, error;

  // Calculates initial conditions based on Taylor expansion
  double aux = 27*u*pow(K , 3)*pow(a , -4);
  xk = pow(aux , 0.25);
  if (xk < 0.5) xk = pow(aux , 1.0/3.0);

  error = 1000;
  // Newton-Raphson algorithm
  while ( error > tol) {
      f = pow(xk , 3.0)*(xk + 1.0) - aux;
      fprime = 3*pow(xk , 2.0) + 4*pow(xk , 3.0);

      xkplusone = xk - f/fprime;

      error = abs(xk-xkplusone);
      xk = xkplusone;
  }

  return xkplusone;

}

} // namespace FLRWBackground
