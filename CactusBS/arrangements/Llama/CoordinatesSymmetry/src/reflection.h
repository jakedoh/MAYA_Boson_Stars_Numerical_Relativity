#ifndef COORDINATESSYMMETRY_H
#define COORDINATESSYMMETRY_H

#include "cctk.h"
#include "cctk_Arguments.h"

void
CoordinatesSymmetry_Apply (CCTK_ARGUMENTS);

void 
CoordinatesSymmetry_GetManualParities(int table, int gi, int *parities);

#endif /* ! defined COORDINATESSYMMETRY_H */
