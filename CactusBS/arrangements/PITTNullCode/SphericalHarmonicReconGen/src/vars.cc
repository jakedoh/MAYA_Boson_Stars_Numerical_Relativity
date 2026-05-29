#include "vars.hh"

namespace SHR {

/* state info vars */
int       initialized = 0;
int       read_data = 0;

int       time_derivative_in_file = 0;

int       starting_iteration = 0;
int       iteration_in_file = 0;
CCTK_REAL initial_time_in_file = 0;

vector<SPH_database*> db(NUM_METRIC_COMPONENTS*3, NULL);

vector<spherical_decomposed_variable<CCTK_REAL>* > C(NUM_METRIC_COMPONENTS, NULL);
vector<spherical_decomposed_variable<CCTK_REAL>* > Cr(NUM_METRIC_COMPONENTS, NULL);
vector<spherical_decomposed_variable<CCTK_REAL>* > Ct(NUM_METRIC_COMPONENTS, NULL);

bool use_Condon_Shortley_phase_factor = false;

MPI_Comm comm_world = MPI_COMM_NULL;

}

