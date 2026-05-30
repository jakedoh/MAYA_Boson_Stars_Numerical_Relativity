
/* Include both levels of constants */
#include "constants.h"
#include "rnsv2/consts.h"

#include "rnsv2/equil.h"
#include "rnsv2/equil_util.h"
#include "rnsv2/nrutil.h"

#ifdef SQR(a)
#undef SQR(a)
#endif

#define SQR(x) ((x)*(x))

#define JUNK 100000

// In rns.c
void RNS_output_state( double r_ratio, double central_energy_density, double mass, double mass0, double Re,
		double Omega, double Omega_K, double J, int with_header);
void RNS_output_polystate( double r_ratio, double central_energy_density, double mass, double mass0, double rns_radius, 
  double Omega, double Omega_K, double J, double eos_Gamma, double eos_Kappa, int with_header);

/********************/
/* In rns_methods.c */
/********************/

/* Create stars with a given angular momentum */ 
void RNS_angmom( double s_gp[SDIV+1], double mu[MDIV+1], double log_e_tab[201], double log_p_tab[201],
  double log_h_tab[201], double log_n0_tab[201], int n_eosentries, char eos_type[],
  double **rho, double **gama, double **alpha, double **omega, double **energy, double **pressure, double **enthalpy, 
  double **velocity_sq, double *v_plus, double *v_minus, double p_center, double h_center, double e_surface, double enthalpy_min,
  double *r_e, double *J, double *final_oblateness, double *final_omega );

/* Create stars for a specified Omega/Omega_K */
void RNS_Omega( double s_gp[SDIV+1], double mu[MDIV+1], double log_e_tab[201], double log_p_tab[201],
  double log_h_tab[201], double log_n0_tab[201], int n_eosentries, char eos_type[],
  double **rho, double **gama, double **alpha, double **omega, double **energy, double **pressure, double **enthalpy, 
  double **velocity_sq, double *v_plus, double *v_minus, double p_center, double h_center, double e_surface, double enthalpy_min,
  double *r_e, double *J, double *final_oblateness, double *final_Omega );

