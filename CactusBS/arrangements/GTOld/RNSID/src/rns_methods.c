// RNSID: rns_methods.c
// Methods for creating a rotating neutron star

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "rns.h"
#define DSIGN(x,y) ( fabs(x-y) )

static inline CCTK_INT mysign( double x) {
  CCTK_REAL realsgn=x/fabs(x);
  return (int) realsgn;
}

/* Create stars of increasing oblateness (by decreasing r_ratio), starting
 * from a non-rotating star, until the angular momentum of the generated star is
 * as specified. */
void RNS_angmom( double s_gp[SDIV+1], double mu[MDIV+1], double log_e_tab[201], double log_p_tab[201],
  double log_h_tab[201], double log_n0_tab[201], int n_eosentries, char eos_type[],
  double **rho, double **gama, double **alpha, double **omega, double **energy, double **pressure, double **enthalpy, 
  double **velocity_sq, double *v_plus, double *v_minus, double p_center, double h_center, double e_surface, double enthalpy_min,
  double *r_e, double *J, double *final_oblateness, double *final_omega )
{

  DECLARE_CCTK_PARAMETERS;
  double rns_radius, Mass, Mass_0, Omega, Omega_K, a_check=0;
  double r_ratio = *final_oblateness;

  // Solver solves in K=1 units. Assume rns_angular_momentum requested in geometric units with K=/=1
  // Base relevant conversion factor is K^(n/2)
  //   J' = K^{-n} J
  double eos_n = 1./(eos_Gamma-1.);
  double kapfactor = pow( eos_Kappa, eos_n );

  double diff_J = rns_angular_momentum/kapfactor - *J;
  double old_diff_J = diff_J;
  if (verbose)
	  CCTK_VInfo(CCTK_THORNSTRING,"Starting iterative loop, diff_J = %g, (J=%g, J_target=%g).", diff_J, *J, rns_angular_momentum);
  while( (diff_J > 0) ) { /* Approach J from below, oblateness from above */

	  /* Find the interval of r_ratio with the correct angular momentum */
	  r_ratio -= delta_oblate;

	  /* Compute the star of this r_ratio */
	  spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, h_center, enthalpy_min, rho, gama, alpha, omega, 
		energy, pressure, enthalpy, velocity_sq, a_check, rns_accuracy, 
		convergence_factor, r_ratio, r_e, &Omega);
	  mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, rho, gama, alpha, omega, energy, pressure, enthalpy, 
		velocity_sq, r_ratio, e_surface, *r_e, Omega, &Mass, &Mass_0, J, &rns_radius, v_plus, 
	        v_minus, &Omega_K);

  	  if ( CCTK_EQUALS(eos_info,"poly")==0)
		RNS_output_state(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, 0);
	  else
		RNS_output_polystate(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, eos_Gamma, eos_Kappa, 0);

	  old_diff_J = diff_J;
	  diff_J = rns_angular_momentum/kapfactor - *J;
	  if ( verbose )
	  	CCTK_VInfo(CCTK_THORNSTRING,"     Computed new J: %g  (old_diff_J = %g, new diff_J = %g).",*J,old_diff_J, diff_J);

  }

  /* The desired star lies between r_ratio and r_ratio + delta_oblateness */
  if (verbose)
	CCTK_VInfo(CCTK_THORNSTRING,"Desired star lies between oblateness of %g and %g. Zooming in.",r_ratio, r_ratio+delta_oblate);

  double oblow, obhigh, flow, fhigh, obmed, fmed;
  oblow = r_ratio;
  obhigh = r_ratio + delta_oblate;
  flow = diff_J;
  fhigh = old_diff_J;

  /* Use Ridder's method to find the correct star (Taken from NR: rework!) */
  /*   f(x) = dJ(oblateness)                                               */
  double cur_oblate, obnew, fnew, sroot;
  if (( flow>0.0 && fhigh<0.0 ) || (flow<0.0 && fhigh>0.0) ) {
     cur_oblate=-1.11e30;
     for ( int iter=1;iter<=max_obits;iter++) {
	  obmed=0.5*(oblow+obhigh);
	  r_ratio = obmed;

	  spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, h_center, enthalpy_min, rho, gama, alpha, omega, 
		energy, pressure, enthalpy, velocity_sq, a_check, rns_accuracy, 
		convergence_factor, r_ratio, r_e, &Omega);
	  mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, rho, gama, alpha, omega, energy, pressure, enthalpy, 
		velocity_sq, r_ratio, e_surface, *r_e, Omega, &Mass, &Mass_0, J, &rns_radius, v_plus, 
	        v_minus, &Omega_K);

	  fmed = -(*J - rns_angular_momentum/kapfactor);
	  sroot=sqrt( SQR(fmed) - flow*fhigh);
	  if (sroot == 0.0) {
		  r_ratio = cur_oblate;
		  break;
	  }

	  obnew = obmed + (obmed-oblow)*( (flow >= fhigh ? 1.0 : -1.0)*fmed/sroot);
	  if (fabs(obnew-cur_oblate) <= radial_accuracy) {
		  r_ratio = cur_oblate;
		  break;
	  }
	  cur_oblate = obnew;
	  r_ratio = cur_oblate;

	  spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, h_center, enthalpy_min, rho, gama, alpha, omega, 
		energy, pressure, enthalpy, velocity_sq, a_check, rns_accuracy, 
		convergence_factor, r_ratio, r_e, &Omega);
	  mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, rho, gama, alpha, omega, energy, pressure, enthalpy, 
		velocity_sq, r_ratio, e_surface, *r_e, Omega, &Mass, &Mass_0, J, &rns_radius, v_plus, 
	        v_minus, &Omega_K);

	  if ( CCTK_EQUALS(eos_info,"poly")==0)
		RNS_output_state(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, 0);
	  else
		RNS_output_polystate(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, eos_Gamma, eos_Kappa, 0);

	  fnew = rns_angular_momentum/kapfactor - *J;
	  if (fnew == 0.0){
		r_ratio = cur_oblate;
		break;
	  }

          if (verbose)
             CCTK_VInfo(CCTK_THORNSTRING,"   f(%g)=%g  f(%g)=%g  f(%g)=%g ... f(%g)=%g", oblow, flow, obmed, fmed, obhigh, fhigh, obnew, fnew );

          // Is root above or below cur_oblate?
          if ( fabs(fmed-fnew) > fabs(fmed) ) { // Root between obnew and obmed

             if ( obnew > obmed ) {
                oblow=obmed;
                flow=fmed;
                obhigh=cur_oblate;
                fhigh=fnew;
             } else {
                oblow=cur_oblate;
                flow=fnew;
                obhigh=obmed;
                fhigh=fmed;
             }

          } else if ( fabs(fhigh-fnew) > fabs(fhigh) ) { // Root is closer to fhigh from fnew
                oblow=cur_oblate;
                flow=fnew;
          } else if ( fabs(flow-fnew) > fabs(flow) ) { // Root is closer to flow from fnew
                obhigh=cur_oblate;
                fhigh=fnew;
          } else if ( fnew == 0 ) {// Somehow landed right on the button? 
                r_ratio = cur_oblate;
                break;
          } else nrerror("Never get here.");

  	  if ( fabs(obhigh-oblow) <= radial_accuracy) { // Accurate enough in radius.
		r_ratio = cur_oblate;
		break;
	  }
     }
  } else {
	  if (fhigh == 0.0)
		  r_ratio += delta_oblate;
	  nrerror("Root must be bracketed.");
  }

  ( *final_oblateness ) = r_ratio;
  ( *final_omega ) = Omega;

}

/* Create stars of increasing oblateness (by decreasing r_ratio), starting
 * from a non-rotating star, until the angular momentum of the generated star is
 * as specified. */
void RNS_Omega( double s_gp[SDIV+1], double mu[MDIV+1], double log_e_tab[201], double log_p_tab[201],
  double log_h_tab[201], double log_n0_tab[201], int n_eosentries, char eos_type[],
  double **rho, double **gama, double **alpha, double **omega, double **energy, double **pressure, double **enthalpy, 
  double **velocity_sq, double *v_plus, double *v_minus, double p_center, double h_center, double e_surface, double enthalpy_min,
  double *r_e, double *J, double *final_oblateness, double *final_omega )
{

  DECLARE_CCTK_PARAMETERS;
  double rns_radius, Mass, Mass_0, Omega, Omega_K, a_check=0;
  double r_ratio = *final_oblateness;
	
  double diff_fOmega = rns_omega_frac_of_Kepler; // We're starting with a non-rotatin star.
  double old_diff_fOmega = diff_fOmega;
  if (verbose)
	  CCTK_VInfo(CCTK_THORNSTRING,"Starting iterative loop, diff_fOmega = %g, (J=%g, Omega_target = %g*Omega_K).", diff_fOmega, *J, rns_omega_frac_of_Kepler);
  while( (diff_fOmega > 0) ) { /* Approach J from below, oblateness from above */

	  /* Find the interval of r_ratio with the correct angular momentum */
	  r_ratio -= delta_oblate;

	  /* Compute the star of this r_ratio */
	  spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, h_center, enthalpy_min, rho, gama, alpha, omega, 
		energy, pressure, enthalpy, velocity_sq, a_check, rns_accuracy, 
		convergence_factor, r_ratio, r_e, &Omega);
	  mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, rho, gama, alpha, omega, energy, pressure, enthalpy, 
		velocity_sq, r_ratio, e_surface, *r_e, Omega, &Mass, &Mass_0, J, &rns_radius, v_plus, 
	        v_minus, &Omega_K);

  	  if ( CCTK_EQUALS(eos_info,"poly")==0)
		RNS_output_state(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, 0);
	  else
		RNS_output_polystate(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, eos_Gamma, eos_Kappa, 0);

	  old_diff_fOmega = diff_fOmega;
	  diff_fOmega = rns_omega_frac_of_Kepler*Omega_K - Omega;
	  if ( verbose )
	  	CCTK_VInfo(CCTK_THORNSTRING,"     Computed new Omega: %g*Omega_K  (old_diff_fOmega = %g, new diff_fOmega = %g).",Omega/Omega_K,old_diff_fOmega, diff_fOmega);

  }

  /* The desired star lies between r_ratio and r_ratio + delta_oblateness */
  if (verbose)
	CCTK_VInfo(CCTK_THORNSTRING,"Desired star lies between oblateness of %g and %g. Zooming in.",r_ratio, r_ratio+delta_oblate);

  double oblow, obhigh, flow, fhigh, obmed, fmed;
  oblow = r_ratio;
  obhigh = r_ratio + delta_oblate;
  flow = diff_fOmega;
  fhigh = old_diff_fOmega;

  /* Use Ridder's method to find the correct star (Taken from NR: rework!) */
  /*   f(x) = dJ(oblateness)                                               */
  double cur_oblate, obnew, fnew, sroot;
  if (( flow>0.0 && fhigh<0.0 ) || (flow<0.0 && fhigh>0.0) ) {
     cur_oblate=-1.11e30;
     for ( int iter=1;iter<=max_obits;iter++) {
	  obmed=0.5*(oblow+obhigh);
	  r_ratio = obmed;

	  spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, h_center, enthalpy_min, rho, gama, alpha, omega, 
		energy, pressure, enthalpy, velocity_sq, a_check, rns_accuracy, 
		convergence_factor, r_ratio, r_e, &Omega);
	  mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, rho, gama, alpha, omega, energy, pressure, enthalpy, 
		velocity_sq, r_ratio, e_surface, *r_e, Omega, &Mass, &Mass_0, J, &rns_radius, v_plus, 
	        v_minus, &Omega_K);

	  fmed = rns_omega_frac_of_Kepler*Omega_K - Omega;
	  sroot=sqrt( SQR(fmed) - flow*fhigh);
	  if (sroot == 0.0) {
		  r_ratio = cur_oblate;
		  break;
	  }

	  obnew = obmed + (obmed-oblow)*( (flow >= fhigh ? 1.0 : -1.0)*fmed/sroot);
	  if (fabs(obnew-cur_oblate) <= radial_accuracy) {
		  r_ratio = cur_oblate;
		  break;
	  }
	  cur_oblate = obnew;
	  r_ratio = cur_oblate;

	  spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, h_center, enthalpy_min, rho, gama, alpha, omega, 
		energy, pressure, enthalpy, velocity_sq, a_check, rns_accuracy, 
		convergence_factor, r_ratio, r_e, &Omega);
	  mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		eos_type, eos_Gamma, rho, gama, alpha, omega, energy, pressure, enthalpy, 
		velocity_sq, r_ratio, e_surface, *r_e, Omega, &Mass, &Mass_0, J, &rns_radius, v_plus, 
	        v_minus, &Omega_K);

	  if ( CCTK_EQUALS(eos_info,"poly")==0)
		RNS_output_state(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, 0);
	  else
		RNS_output_polystate(r_ratio, central_energy_density, Mass, Mass_0, rns_radius, Omega, Omega_K, *J, eos_Gamma, eos_Kappa, 0);

	  fnew = rns_omega_frac_of_Kepler*Omega_K - Omega;
	  if (fnew == 0.0){
		r_ratio = cur_oblate;
		break;
	  }

          if (verbose)
             CCTK_VInfo(CCTK_THORNSTRING,"   f(%g)=%g  f(%g)=%g  f(%g)=%g ... f(%g)=%g", oblow, flow, obmed, fmed, obhigh, fhigh, obnew, fnew );

          // Is root above or below cur_oblate?

          if ( SIGN(fmed,fnew) != fmed ) { 
             oblow=obmed;
             flow=fmed;
             obhigh=cur_oblate;
             fhigh=fnew;
          } else if (SIGN(flow,fnew) != flow ) {
             obhigh=cur_oblate;
             fhigh=fnew;
          } else if (SIGN(fhigh,fnew) != fhigh ) {
             oblow=cur_oblate;
             flow=fnew;
          } else nrerror("never get here.");
          if ( fabs(obhigh-oblow) <= radial_accuracy )  {
             r_ratio = cur_oblate;
             break;
          }

     }

  } else {
	  if (fhigh == 0.0)
		  r_ratio += delta_oblate;
	  nrerror("Root must be bracketed.");
  }

  ( *final_oblateness ) = r_ratio;
  ( *final_omega ) = Omega;

}

