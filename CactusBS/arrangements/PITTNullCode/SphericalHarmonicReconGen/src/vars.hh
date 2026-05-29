#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <math.h>
#include <complex>

#ifndef _VARS__HHHH_
#define _VARS__HHHH_

#include "sph_database.hh"

#define NUM_METRIC_COMPONENTS 10


namespace SHR {

/* state info vars */
extern int       initialized;
extern int       read_data;

extern int       time_derivative_in_file;

extern int       starting_iteration;
extern int       iteration_in_file;
extern CCTK_REAL initial_time_in_file;

extern vector<SPH_database*> db; 

extern vector<spherical_decomposed_variable<CCTK_REAL>* > C;
extern vector<spherical_decomposed_variable<CCTK_REAL>* > Cr;
extern vector<spherical_decomposed_variable<CCTK_REAL>* > Ct;

extern bool use_Condon_Shortley_phase_factor;

extern MPI_Comm comm_world;

}


#endif