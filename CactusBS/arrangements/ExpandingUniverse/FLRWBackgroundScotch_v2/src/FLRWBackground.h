#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <util_Table.h>
#include <Symmetry.h>

#include "Kranc.hh"

#define PI 3.14159265358979323846264338328
#define TINY 0.0000000000000000001

extern "C" void initialHydro_FLRW(CCTK_ARGUMENTS);

extern "C" void initialMetric_FLRW(CCTK_ARGUMENTS);

extern "C" void FLRWBackground_SelectBoundConds(CCTK_ARGUMENTS);

extern "C" void FLRWBackground_HydroBoundary(CCTK_ARGUMENTS);

extern "C" void FLRWBackground_HydroBoundaryConformal(CCTK_ARGUMENTS);

extern "C" void FLRWBackground_MetricBoundary(CCTK_ARGUMENTS);

extern "C" void FLRWBackground_MetricBoundaryConformal(CCTK_ARGUMENTS);

extern "C" void FLRWBackground_RegisterSymmetries(CCTK_ARGUMENTS);

extern "C" double getEps(double u , double K , double a, double tol);
