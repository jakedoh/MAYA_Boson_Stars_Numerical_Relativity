// RNSID: rns.c
// Main code for interfacing with Stergioulas' RNS code from Cactus.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "cctk_Loop.h"
#include "rns.h"
#include "HydroBase.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#define INITVAL 1000

#define DIM 2
#define INTERP_NPOINTS 1
#define NUM_INPUT_ARRAYS 8
#define BASE_OUTPUT_ARRAYS 8
#define NUM_TO_D 4
#define NUM_OUTPUT_ARRAYS 8+DIM*NUM_TO_D
#define CRMIN 1e-3

#define SURF_START_PRESS_F  1000
#define SURF_INTERP_NPOINTS 100
#define SURF_NUM_INPUT_ARRAYS 2
#define SURF_NUM_OUTPUT_ARRAYS 2

#define IDX2D(i_,j_) ( (i_) + npts_intrad*(j_) ) 

/* Local function declarations */
static CCTK_INT RNS_output_solution( CCTK_REAL *rho, CCTK_REAL *gama, CCTK_REAL *alpha, CCTK_REAL *omega, CCTK_REAL *energy, 
           CCTK_REAL *pressure, CCTK_REAL *enthalpy, CCTK_REAL *velocity_sq, 
           CCTK_REAL ns_radius, CCTK_REAL ns_radius_iso, CCTK_REAL ns_mass, CCTK_REAL ns_oblateness, CCTK_REAL ns_omega );
static CCTK_INT RNS_input_solution( CCTK_REAL *rho, CCTK_REAL *gama, CCTK_REAL *alpha, CCTK_REAL *omega, CCTK_REAL *energy, 
           CCTK_REAL *pressure, CCTK_REAL *enthalpy, CCTK_REAL *velocity_sq, CCTK_REAL *ns_radius, CCTK_REAL *ns_radius_iso, 
           CCTK_REAL *ns_mass, CCTK_REAL *ns_oblateness, CCTK_REAL *ns_omega );
static CCTK_INT perturb_star( const CCTK_INT pert_idx, const CCTK_REAL r_rs, const CCTK_REAL mu, CCTK_REAL xx, CCTK_REAL yy, CCTK_REAL r_inv,  CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *press );
static CCTK_REAL *rns_allocate_array(CCTK_INT npoints, const char *name);
static void fill_GA_ghostzones( CCTK_REAL * rnsid_rho, CCTK_REAL * rnsid_gama, CCTK_REAL * rnsid_alpha, CCTK_REAL * rnsid_omega, 
           CCTK_REAL * rnsid_energy, CCTK_REAL * rnsid_pressure, CCTK_REAL * rnsid_enthalpy, 
           CCTK_REAL * rnsid_velocity_sq );

/* External functions */
void RNS_InitGA(CCTK_ARGUMENTS);
void Solve_RNS(CCTK_ARGUMENTS);
void Set_RNS_ID(CCTK_ARGUMENTS);

/* Global variable */
static int firstcall=1;

/******************/
/* Main functions */
/******************/
void Solve_RNS(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* We will eventually copy the solution of the 8 quantities into larger 
     arrays with NGHOSTS ghost zones filled with the appropriate symmetries */
  int npts_intrad = radial_spokes+2*nghosts;
  int npts_intang = angular_spokes+2*nghosts;

  /* Setup EOS and defaults */
  double log_e_tab[201], log_p_tab[201], log_h_tab[201], log_n0_tab[201];
  int n_eosentries;
  double enthalpy_min, e_surface, p_surface;
  char eos_type[5];
  double eos_n = 1./(eos_Gamma-1.);
  double kapn = 1.; 
  double ikapn = 1.;

  // Solver solves in K=1 units. Output in K =/= 1 units
  // Base relevant conversion factor is K^n


  if ( CCTK_EQUALS(eos_info,"poly") == 1) {
     sprintf(eos_type,"poly");
     e_surface = 0.;
     p_surface = 0.;
     enthalpy_min = 0.;
     kapn = pow( eos_Kappa, eos_n );
     ikapn = 1./kapn;
  } else {
     CCTK_WARN(0,"You specify a non-polytype NS. While we can solve for this, "
                 "the filling of the hydro variables currently assumes polytype." );
     sprintf(eos_type,"tab");
     load_eos( eos_type, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_eosentries);
     e_surface = 7.8*SQR(C)*KSCALE;
     p_surface = 1.01e8*KSCALE;
     enthalpy_min = 1.0/(SQR(C));
  }

  if ( *ns_mass > 0 ) {
     firstcall=0;
  }

  /* Should only solve on head node */
  //if ( !read_rns_from_file || ( firstcall>0 && CCTK_MyProc(cctkGH)==0) )
  if ( read_rns_from_file==0 && (firstcall>0) )
  {
     /* Set up grid for solving (not main Cactus grid!) */
     double mu[angular_spokes+1], s_gp[radial_spokes+1]; 
     make_grid( s_gp, mu );

     double **rns_rho, **rns_gama, **rns_omega, **rns_alpha;
     rns_rho = dmatrix(1, radial_spokes, 1, angular_spokes);
     rns_gama = dmatrix(1, radial_spokes, 1, angular_spokes);
     rns_alpha = dmatrix(1, radial_spokes, 1, angular_spokes);
     rns_omega = dmatrix(1, radial_spokes, 1, angular_spokes);

     double **rns_energy, **rns_pressure, **rns_enthalpy, **rns_velocity_sq;
     rns_energy = dmatrix(1, radial_spokes, 1, angular_spokes);
     rns_pressure = dmatrix(1, radial_spokes, 1, angular_spokes);
     rns_enthalpy = dmatrix(1, radial_spokes, 1, angular_spokes);
     rns_velocity_sq = dmatrix(1, radial_spokes, 1, angular_spokes);

     /* We keep a local copy here, then assign to grid arrays later */
     double *v_plus, *v_minus;
     v_plus = dvector(1, radial_spokes);
     v_minus = dvector(1, radial_spokes);

     /* Calculate the pressure and enthalpy at the center of the star.
      * This will be asserted upon iterations. */
     double p_center, h_center;
     make_center( eos_info, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
	eos_type, eos_Gamma, kapn*central_energy_density, &p_center, &h_center );
     if (verbose)
	  CCTK_VInfo(CCTK_THORNSTRING,"Set up the central NS measures: (e, P, h) = (%g,%g,%g), or for K=/=1 polytrope (e,P,h) = (%g,%g,%g)", 
                     kapn*central_energy_density, p_center, h_center, central_energy_density, ikapn*p_center, h_center );

     /* Compute a spherical star as a first guess for the rotating star */
     double Omega, r_e; /* r_e is used to convert between s and r */
     sphere( s_gp, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
	eos_type, eos_Gamma, kapn*central_energy_density, p_center, h_center, p_surface, e_surface,
	rns_rho, rns_gama, rns_alpha, rns_omega, &r_e);
     //if (verbose)
     //  CCTK_VInfo(CCTK_THORNSTRING,"Found a non-rotating star with these central quantities whose radius is %g (Schw).", r_e);

     /* Compute the metric of a star with a given oblateness, where oblateness is 
      * given by an index r_ratio of the length of the axis connecting the centre
      * of the star to one of the poles to the radius of the star's equator. 
      *
      * The metric is then given through the functions alpha, rho, gama, and omega
      * as given in the documentation. */
     double a_check=0;
     double r_ratio = 1.;
     spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
   	eos_type, eos_Gamma, h_center, enthalpy_min, rns_rho, rns_gama, rns_alpha, rns_omega,
	rns_energy, rns_pressure, rns_enthalpy, rns_velocity_sq, a_check, rns_accuracy, convergence_factor,
	r_ratio, &r_e, &Omega);
     if (verbose)
	  CCTK_VInfo(CCTK_THORNSTRING,"Computed a spinning NS of radius %g, oblateness %g, and angular velocity %g.", 
                     sqrt(kapn)*r_e, r_ratio, sqrt(ikapn)*Omega);

     /* Compute the various equilibrium quantities */
     double Mass_0, J, Omega_K;
     mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
	eos_type, eos_Gamma, rns_rho, rns_gama, rns_alpha, rns_omega, rns_energy, rns_pressure, rns_enthalpy, 
	rns_velocity_sq, r_ratio, e_surface, r_e, Omega, ns_mass, &Mass_0, &J, ns_radius, v_plus, 
        v_minus, &Omega_K);

     /* Print out the information about the stellar model (non-spinning) */
     if ( CCTK_EQUALS(eos_info,"poly")==0)
	RNS_output_state(r_ratio, central_energy_density, *ns_mass, Mass_0, *ns_radius, Omega, Omega_K, J, 1);
     else
	RNS_output_polystate(r_ratio, central_energy_density, *ns_mass, Mass_0, *ns_radius, Omega, Omega_K, J, 
                             eos_Gamma, eos_Kappa, 1);

     /* Iterate to desired physical star depending on parameter choice ns_search_type. */
     if ( CCTK_EQUALS(ns_search_type,"AngMom") && (rns_angular_momentum > 0) ) {
	RNS_angmom( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		     eos_type, rns_rho, rns_gama, rns_alpha, rns_omega, rns_energy, 
		     rns_pressure, rns_enthalpy, rns_velocity_sq, v_plus, v_minus, p_center,
		     h_center, e_surface, enthalpy_min, &r_e, &J, &r_ratio, &Omega );
     }
     if ( CCTK_EQUALS(ns_search_type,"Omega") && (rns_omega_frac_of_Kepler > 0) ) {
	RNS_Omega( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
		     eos_type, rns_rho, rns_gama, rns_alpha, rns_omega, rns_energy, 
		     rns_pressure, rns_enthalpy, rns_velocity_sq, v_plus, v_minus, p_center,
		     h_center, e_surface, enthalpy_min, &r_e, &J, &r_ratio, &Omega );
     }

     /* Recompute and spit out equilibrium quantities for final star. */
     mass_radius( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, n_eosentries,
	eos_type, eos_Gamma, rns_rho, rns_gama, rns_alpha, rns_omega, rns_energy, rns_pressure, rns_enthalpy, 
	rns_velocity_sq, r_ratio, e_surface, r_e, Omega, ns_mass, &Mass_0, &J, ns_radius, v_plus, 
        v_minus, &Omega_K);


     if ( CCTK_EQUALS(eos_info,"poly")==0)
	RNS_output_state(r_ratio, central_energy_density, *ns_mass, Mass_0, *ns_radius, Omega, Omega_K, J, 0);
     else
	RNS_output_polystate(r_ratio, central_energy_density, *ns_mass, Mass_0, *ns_radius, Omega, Omega_K, J, 
                             eos_Gamma, eos_Kappa, 0);

     /* Set iso version of ns_radius and convert scalars to K =/= 1 units if necessary (kapn=1 otherwise) */
     *ns_radius *= sqrt(kapn);
     *ns_mass *= sqrt(kapn);
     //*ns_radius_iso = 0.5*( *ns_radius - *ns_mass + sqrt( *ns_radius*(*ns_radius-2.*(*ns_mass)) ) );
     *ns_radius_iso = sqrt(kapn)*r_e;
     *ns_oblateness = r_ratio;
     *ns_omega = sqrt(ikapn)*Omega;

     /* Fill in what's already generated, converting to K =/= 1 units */
     for ( int i=1; i<=radial_spokes; i++ ) {
         for ( int j=1; j<=angular_spokes; j++ ){
             CCTK_INT angIdx=IDX2D( (i-1)+nghosts,(j-1)+nghosts); // rns idx starts at 1
   	     rnsid_rho[angIdx]  	 	= rns_rho[i][j];
             rnsid_gama[angIdx] 	 	= rns_gama[i][j];
             rnsid_alpha[angIdx] 	 	= rns_alpha[i][j];
             rnsid_omega[angIdx] 	 	= sqrt(ikapn)*rns_omega[i][j];
             rnsid_energy[angIdx]  	        = ikapn*rns_energy[i][j];
             rnsid_pressure[angIdx] 	        = ikapn*rns_pressure[i][j];
             rnsid_enthalpy[angIdx]   	        = rns_enthalpy[i][j];
             rnsid_velocity_sq[angIdx]  	= rns_velocity_sq[i][j];
         }
     }
     free_dmatrix(rns_rho,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_gama,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_omega,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_alpha,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_energy,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_enthalpy,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_pressure,1,radial_spokes,1,angular_spokes); 
     free_dmatrix(rns_velocity_sq,1,radial_spokes,1,angular_spokes); 

     fill_GA_ghostzones( rnsid_rho, rnsid_gama, rnsid_alpha, rnsid_omega, 
        rnsid_energy, rnsid_pressure, rnsid_enthalpy, rnsid_velocity_sq );

     /* Write out RNS solution for debugging or future use 		* 
      *    Rest of the processors and go on and start filling in the grid 	*/
     if ( write_out_rns && CCTK_MyProc(cctkGH) == 0 ) {
 	int write_err = RNS_output_solution( rnsid_rho, rnsid_gama, rnsid_alpha, rnsid_omega, rnsid_energy, 
                               rnsid_pressure, rnsid_enthalpy, rnsid_velocity_sq, *ns_radius, *ns_radius_iso, 
                               *ns_mass, *ns_oblateness, *ns_omega);
	if ( write_err < 0 )
		CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Error in writing solution: Status %d received.", write_err);
     }


  }
  if ( read_rns_from_file && firstcall>0 ) {

     /* If we didn't create the RNS this simulation (or not on this process), we 
        need to read it in now.  Read data arrays from file (with ghosts zones)  */
     int input_err = RNS_input_solution( rnsid_rho, rnsid_gama, rnsid_alpha, rnsid_omega,
		      rnsid_energy, rnsid_pressure, rnsid_enthalpy, rnsid_velocity_sq, ns_radius, ns_radius_iso, 
                      ns_mass, ns_oblateness, ns_omega );
     if ( input_err < 0 )
        CCTK_WARN(0,"RNS Reading error.  Should never get here (should have died within the function).");

     /* Write out solution, just to verify read-in */
     if ( write_out_rns && CCTK_MyProc(cctkGH) == 0 ) {
       int write_err = RNS_output_solution( rnsid_rho, rnsid_gama, rnsid_alpha, rnsid_omega, rnsid_energy, 
                               rnsid_pressure, rnsid_enthalpy, rnsid_velocity_sq, *ns_radius, *ns_radius_iso, 
                               *ns_mass, *ns_oblateness, *ns_omega);
       if ( write_err < 0 )
               CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Error in writing solution: Status %d received.", write_err);
     }

  }

  // Make sure *something* has been set
  assert( *ns_mass>0 );

}


/************/
/* Fill GFs */
/************/
void Set_RNS_ID(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* We will eventually copy the solution of the 8 quantities into larger 
     arrays with NGHOSTS ghost zones filled with the appropriate symmetries */
  int npts_intrad = radial_spokes+2*nghosts;
  int npts_intang = angular_spokes+2*nghosts;

  /* Setup EOS and defaults */
  double log_e_tab[201], log_p_tab[201], log_h_tab[201], log_n0_tab[201];
  int n_eosentries;
  double enthalpy_min, e_surface, p_surface;
  char eos_type[5];

  if ( CCTK_EQUALS(eos_info,"poly") == 1) {
     sprintf(eos_type,"poly");
     e_surface = 0.;
     p_surface = 0.;
     enthalpy_min = 0.;
  } else {
     CCTK_WARN(0,"You specify a non-polytype NS. While we can solve for this, "
                 "the filling of the hydro variables currently assumes polytype." );
     sprintf(eos_type,"tab");
     load_eos( eos_type, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_eosentries);
     e_surface = 7.8*SQR(C)*KSCALE;
     p_surface = 1.01e8*KSCALE;
     enthalpy_min = 1.0/(SQR(C));
  }

  CCTK_INT set_lapse=0;
  CCTK_INT set_shift=0;
  if ( CCTK_EQUALS(initial_lapse,"RNS") )
     set_lapse=1;
  if ( CCTK_EQUALS(initial_shift,"RNS") )
     set_shift=1;

  if (verbose)
     CCTK_VInfo(CCTK_THORNSTRING,"Setting up GFs given RNS Solution.");

  /* Set perturbation info */
  CCTK_INT pert_idx = -1;
  if ( add_pert ) {

     if ( CCTK_EQUALS(pert_type,"Radial Zink") ){
        pert_idx = 1;
     } 
     if ( CCTK_EQUALS(pert_type,"f2 Zink") ) {
        pert_idx = 2;
     }
  }

  /* Set atmo info */
  CCTK_REAL pmin, epsmin;
  pmin = eos_Kappa*pow(rho_atmosphere,eos_Gamma);
  epsmin = pmin/( rho_atmosphere*(eos_Gamma-1) ); 


  /*****************************************************************************/
  /** This function is called once per level, so interpolate once per level   **/
  /* 1. Set up Handles, including setting up the 2D coordinate system          */
  /* 2. Create input arrays of (s,mu) coordinates given level's x,y,z data     */
  /* 3. Interpolate rnsid_* onto temporary arrays, all at once                 */
  /* 4. With OMP, loop over locations and populate gridfunctions from rnsid_*  */
  /*****************************************************************************/

  /* Handles for interpolation onto grid */
  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);
  const int param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if (param_table_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  }

  #define DERIV(x) x
  CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS];
  CCTK_INT opcodes[NUM_OUTPUT_ARRAYS];
  for(int idx = 0 ; idx < BASE_OUTPUT_ARRAYS ; idx++) {
    operand_indices[idx] = idx;
    opcodes[idx] = DERIV(0);
  }
  for(int idx = BASE_OUTPUT_ARRAYS; idx < BASE_OUTPUT_ARRAYS+NUM_TO_D; idx++) {
    operand_indices[idx] = idx - BASE_OUTPUT_ARRAYS; // Operate on metric potentials
    opcodes[idx] = DERIV(2); // df/dx^(1)
  }
  for(int idx = BASE_OUTPUT_ARRAYS+NUM_TO_D; idx < BASE_OUTPUT_ARRAYS+2*NUM_TO_D; idx++) {
    operand_indices[idx] = idx - (BASE_OUTPUT_ARRAYS+NUM_TO_D);  // Operate on metric potentials
    opcodes[idx] = DERIV(1); // df/dx^(2)
  }
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operand_indices, "operand_indices");
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        opcodes, "operation_codes");


  /* Input coordinates of level gridpoints */
  CCTK_INT npts_lvl, nx, ny, nz;
  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];
  npts_lvl = nx*ny*nz; 

  /* Specify resolution in (s,mu) space */
  const CCTK_REAL ds = SMAX / (SDIV - 1.0);
  const CCTK_REAL dmu =  1. / (MDIV - 1.0);

  CCTK_REAL * rns_alpha = rns_allocate_array( npts_lvl, "rns_alpha" ); 
  CCTK_REAL * rns_rho   = rns_allocate_array( npts_lvl, "rns_rho" ); 
  CCTK_REAL * rns_gama  = rns_allocate_array( npts_lvl, "rns_gama" ); 
  CCTK_REAL * rns_omega = rns_allocate_array( npts_lvl, "rns_omega" ); 
  CCTK_REAL * rns_e     = rns_allocate_array( npts_lvl, "rns_e" ); 
  CCTK_REAL * rns_press = rns_allocate_array( npts_lvl, "rns_press" ); 
  CCTK_REAL * rns_h     = rns_allocate_array( npts_lvl, "rns_h" ); 
  CCTK_REAL * rns_v2    = rns_allocate_array( npts_lvl, "rns_v2" ); 

  CCTK_REAL * dalpha_dmu = rns_allocate_array( npts_lvl, "dalpha_dmu" ); 
  CCTK_REAL * drho_dmu   = rns_allocate_array( npts_lvl, "drho_dmu" ); 
  CCTK_REAL * dgama_dmu  = rns_allocate_array( npts_lvl, "dgama_dmu" ); 
  CCTK_REAL * domega_dmu = rns_allocate_array( npts_lvl, "domega_dmu" ); 

  CCTK_REAL * dalpha_ds = rns_allocate_array( npts_lvl, "dalpha_ds" ); 
  CCTK_REAL * drho_ds   = rns_allocate_array( npts_lvl, "drho_ds" ); 
  CCTK_REAL * dgama_ds  = rns_allocate_array( npts_lvl, "dgama_ds" ); 
  CCTK_REAL * domega_ds = rns_allocate_array( npts_lvl, "domega_ds" ); 

  CCTK_REAL * s_coords  = rns_allocate_array( npts_lvl, "s_coords" ); 
  CCTK_REAL * mu_coords = rns_allocate_array( npts_lvl, "mu_coords" ); 


  //#pragma omp parallel for 
  //for ( int k=0; k<nz; k++ )
  //  for ( int j=0; j<ny; j++ )
  //   for ( int i=0; i<nx; i++ ) {
  CCTK_LOOP3_ALL(rnsid_coordsetup,cctkGH,i,j,k) {

      CCTK_INT in_idx = i + nx*(j+k*ny);
      CCTK_INT gf_idx = CCTK_GFINDEX3D( cctkGH, i,j,k );

      CCTK_REAL zL, rL,xL,yL;
      zL = z[gf_idx];
      xL = x[gf_idx];
      yL = y[gf_idx];
       
      rL = sqrt( SQR(x[gf_idx]) + SQR(y[gf_idx]) + SQR(zL) );
      
      CCTK_REAL ctheta = ( rL > CRMIN ) ? zL/rL : 0;
      CCTK_REAL ctheta2 = (rL > CRMIN ) ? 1./rL : 0;
      CCTK_REAL sL = rL / (rL + (*ns_radius_iso) );
      CCTK_REAL muL = sqrt(SQR(ctheta)); // Solver coordinates range mu in [0:1]
      CCTK_REAL muL2 = sqrt(SQR(ctheta2)); 
      s_coords[in_idx]  = sL;    
      mu_coords[in_idx] = muL;    

  }
  CCTK_ENDLOOP3_ALL(rnsid_coordsetup);

  /* NOTE: Interpolator coordinates switched: Interpolator expecting Fortran ordering when array allocated with C ordering */
  const CCTK_REAL coord_origin[2] = { -ds*nghosts, -dmu*nghosts };
  const CCTK_REAL coord_delta[2] = { ds, dmu };
  const void *interp_coords[2] = { (const void *) s_coords, (const void *) mu_coords };

  // Input/Output Array Dimension/Type Type info
  const CCTK_INT input_array_dims[DIM] = { npts_intrad, npts_intang };  /* Interpolator Coordinates Switched for ordering issue */
  CCTK_INT input_array_types[NUM_INPUT_ARRAYS];
  for(int i = 0 ; i < NUM_INPUT_ARRAYS ; i++) {
     input_array_types[i] = CCTK_VARIABLE_REAL;
  }

  const void *input_arrays[NUM_INPUT_ARRAYS] = {
     (const void *) rnsid_alpha,
     (const void *) rnsid_rho,
     (const void *) rnsid_gama,
     (const void *) rnsid_omega,
     (const void *) rnsid_energy,
     (const void *) rnsid_pressure,
     (const void *) rnsid_enthalpy,
     (const void *) rnsid_velocity_sq };

  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS ; i++) {
     output_array_types[i] = CCTK_VARIABLE_REAL;
  }


  void * output_arrays[NUM_OUTPUT_ARRAYS]
        = { (void *) rns_alpha,
            (void *) rns_rho,
            (void *) rns_gama,
            (void *) rns_omega,
     
            (void *) rns_e,
            (void *) rns_press,
            (void *) rns_h,
            (void *) rns_v2,
     
            (void *) dalpha_dmu,
            (void *) drho_dmu,
            (void *) dgama_dmu,
            (void *) domega_dmu,
     
            (void *) dalpha_ds,
            (void *) drho_ds,
            (void *) dgama_ds,
            (void *) domega_ds };

  /* Actual interpolation call */
  CCTK_INT ierr = CCTK_InterpLocalUniform(DIM,
                        operator_handle,
                        param_table_handle,
                        coord_origin,
                        coord_delta,
                        npts_lvl,
                          CCTK_VARIABLE_REAL,
                          interp_coords,
                        NUM_INPUT_ARRAYS, // Number of input arrays
                          input_array_dims,
                          input_array_types,
                          input_arrays,
                        NUM_OUTPUT_ARRAYS, // Number of output arrays
                          output_array_types,
                          output_arrays);

  if (ierr != 0) {
     CCTK_WARN(0,"RNS Interpolation onto level grid screwed up.");
  }
  Util_TableDestroy(param_table_handle);


  /* Create a universal r_surface(theta) array for future reference, but only if oblateness != 1 */
  CCTK_REAL rsurf[MDIV];
  CCTK_REAL mu_axis[MDIV];
  if ( *ns_oblateness != 1. ) {

     CCTK_REAL dmu_surf = dmu;

     /* Hand-made shooter method for each mu in the solver grid. */
     /* Inefficient, but not likely to be a time sink            */
     for ( int mu_idx = 0; mu_idx < MDIV; mu_idx++ ) {

         /* Find a place to start, starting at a certain amount above pmin */
         CCTK_REAL p_at_mu[SDIV];
         CCTK_REAL s_at_mu[SDIV];
         for ( int i=0; i<SDIV; i++ ) {
             CCTK_INT idx = IDX2D( i ,mu_idx);
             p_at_mu[i] = rnsid_pressure[idx];
             s_at_mu[i] = i*ds; 
         }
         CCTK_INT low_i = (int)SDIV/4;
         hunt( p_at_mu, SDIV, SURF_START_PRESS_F*pmin, &low_i );

         /* More refined interpolation */
         CCTK_REAL s_surf = -1; 
         CCTK_INT  s_idx = low_i;
         while ( s_idx<SDIV && s_surf<0 ) {
            if ( p_at_mu[s_idx] <= pmin ) {
               s_surf =  s_at_mu[s_idx-1];
            } else {
               s_idx++;
            }
         }

         /* Convert s to r */
         rsurf[mu_idx] = (*ns_radius_iso) * (s_surf/(1.-s_surf)); 
         mu_axis[mu_idx] = dmu*mu_idx;

         /* Check */
         if (s_surf<0 || rsurf[mu_idx] <0.001 ) {
	    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Error finding surface at mu[%d]=%g, s_surf=%g. low_i=%d, p(low_i)=%g, pmin=%g. rsurf=%g",
                       mu_idx, mu_axis[mu_idx], s_surf, low_i, p_at_mu[low_i], pmin, rsurf[mu_idx]);
         }

     }

  } else {
    /* Fill it in anyways */
    for ( int i=0; i<MDIV; i++ ) {
        rsurf[i] = *ns_radius_iso;
    }
  }

  //CCTK_VInfo(CCTK_THORNSTRING,"Here.");
  char pert_fname[50];
  FILE *out_pert;
  if ( add_pert ) { 
     sprintf( pert_fname, "output_perturbation_%d.asc", CCTK_MyProc(cctkGH) );
     out_pert = fopen( pert_fname,"w"); 
  }

  CCTK_LOOP3_ALL(setRNSID, cctkGH, i,j,k ) {

		CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
		CCTK_REAL xL, yL, zL, rL, rbarL;
		xL = x[idx];
		yL = y[idx];
		zL = z[idx];
		rL = sqrt( SQR(xL) + SQR(yL) + SQR(zL) );
		rbarL = sqrt( SQR(xL) + SQR(yL) );

        CCTK_REAL invrbar = ( rbarL > CRMIN ) ? 1./rbarL : 0;
		CCTK_REAL ctheta = ( rL > CRMIN ) ? zL/rL : 0;
		CCTK_REAL ctheta2 = ( rL > CRMIN ) ? 1./rL : 0;
		const CCTK_REAL sL = rL / (rL + *ns_radius_iso);
		const CCTK_REAL muL = sqrt(SQR(ctheta)); // Solver coordinates range mu in [0:1]
		const CCTK_REAL muL2 = sqrt(SQR(ctheta2)); // Solver coordinates range mu in [0:1]

                /* Local rns solution */
                CCTK_INT  int_idx = i + nx*( j + ny*k );
                CCTK_REAL rns_alphaL = rns_alpha[int_idx];
                CCTK_REAL rns_rhoL   = rns_rho[int_idx];
                CCTK_REAL rns_gamaL  = rns_gama[int_idx];
                CCTK_REAL rns_omegaL = rns_omega[int_idx];

                CCTK_REAL rns_energyL   = rns_e[int_idx];
                CCTK_REAL rns_enthalpyL = rns_h[int_idx];
                CCTK_REAL rns_pressureL    = rns_press[int_idx];
                CCTK_REAL rns_velocity_sqL = rns_v2[int_idx];

                CCTK_REAL dalpha_dmuL  = dalpha_dmu[int_idx];
                CCTK_REAL drho_dmuL  = drho_dmu[int_idx];
                CCTK_REAL dgama_dmuL  = dgama_dmu[int_idx];
                CCTK_REAL domega_dmuL  = domega_dmu[int_idx];
     
                CCTK_REAL dalpha_dsL  = dalpha_ds[int_idx];
                CCTK_REAL drho_dsL  = drho_ds[int_idx];
                CCTK_REAL dgama_dsL  = dgama_ds[int_idx];
                CCTK_REAL domega_dsL  = domega_ds[int_idx];
     
		// Fill 3-metric: Standard 3+1 from coordinates in Cartesian. No tricks.
		CCTK_REAL gxxL, gxyL, gxzL, gyyL, gyzL, gzzL;
                CCTK_REAL kxxL, kxyL, kxzL, kyyL, kyzL, kzzL;
                CCTK_REAL alpL, betaxL, betayL, betazL;
                alpL = exp( 0.5*(rns_rhoL + rns_gamaL) );

                CCTK_REAL e2alp, exp_sigma;
		e2alp = exp( 2 * rns_alphaL );
                exp_sigma = exp( rns_gamaL - rns_rhoL );

                // Setup metric.
                if ( rbarL < RMIN ) { // Stationary, simple rotational axis

                   gxxL = e2alp;
                   gyyL = e2alp;
                   gzzL = e2alp;
                   gxzL = 0.;
                   gyzL = 0.;
                   gxyL = 0.;

                   kxxL = 0.;
                   kxyL = 0.;
                   kxzL = 0.;
                   kyyL = 0.;
                   kyzL = 0.;
                   kzzL = 0.;

                   betaxL = 0.;
                   betayL = 0.;
                   betazL = 0.;

                } else {

                   // Let sigma = rns_gama - rns_rho

                   gxxL = ( exp_sigma*SQR(yL) + SQR(xL)*e2alp )/SQR(rbarL);
                   gyyL = ( exp_sigma*SQR(xL) + SQR(yL)*e2alp )/SQR(rbarL); 
                   gzzL = e2alp; 

                   gxyL = ( -exp_sigma*xL*yL + e2alp*xL*yL )/SQR(rbarL);
                   gxzL = 0.;
                   gyzL = 0.;

                   // Coord derivs

                   CCTK_REAL dr_x = xL / rL; 
                   CCTK_REAL dr_y = yL / rL;
                   CCTK_REAL dr_z = zL / rL;

                   CCTK_REAL ds_dx = *ns_radius_iso * dr_x / SQR( rL + *ns_radius_iso );
                   CCTK_REAL ds_dy = *ns_radius_iso * dr_y / SQR( rL + *ns_radius_iso );
                   CCTK_REAL ds_dz = *ns_radius_iso * dr_z / SQR( rL + *ns_radius_iso );

                   CCTK_REAL dmu_dx = -dr_x * zL/SQR(rL); 
                   CCTK_REAL dmu_dy = -dr_y * zL/SQR(rL);
                   CCTK_REAL dmu_dz = ( 1. - SQR(muL) )/rL; 

                   // Derivs of the metric potentials, letting sig=gama-rho 
                   CCTK_REAL domega_dx, domega_dy, domega_dz;
                   domega_dx = domega_dmuL*dmu_dx + domega_dsL*ds_dx; 
                   domega_dy = domega_dmuL*dmu_dy + domega_dsL*ds_dy;
                   if ( zL == 0 ) {
                      domega_dz = 0.;
                   } else {
                      domega_dz = domega_dmuL*dmu_dz + domega_dsL*ds_dz; 
                   }

                   CCTK_REAL dalpha_dx, dalpha_dy, dalpha_dz;
                   dalpha_dx = dalpha_dmuL*dmu_dx + dalpha_dsL*ds_dx;
                   dalpha_dy = dalpha_dmuL*dmu_dy + dalpha_dsL*ds_dy;
                   dalpha_dz = dalpha_dmuL*dmu_dz + dalpha_dsL*ds_dz;

                   CCTK_REAL dsigma_dx, dsigma_dy, dsigma_dz;
                   dsigma_dx = (dgama_dmuL-drho_dmuL)*dmu_dx + (dgama_dsL-drho_dsL)*ds_dx;
                   dsigma_dy = (dgama_dmuL-drho_dmuL)*dmu_dy + (dgama_dsL-drho_dsL)*ds_dy;
                   dsigma_dz = (dgama_dmuL-drho_dmuL)*dmu_dz + (dgama_dsL-drho_dsL)*ds_dz;
                  
                   CCTK_REAL invlapse= exp( -0.5*(rns_rhoL + rns_gamaL) );

                   // See included RNS.nb Mathematica notebook for derivation details
                   //    Simplified to take into account that the solution is axisymmetric,
                   //    meaning kij only depends on derivatives of omega
                   kxxL =      invlapse * yL * exp_sigma * domega_dx
                             + invlapse*SQR(xL/invrbar)*rns_omegaL*e2alp*(yL*dalpha_dx - xL*dalpha_dy) // line vanishes in axisymmetry
                             + invlapse*SQR(yL/invrbar)*rns_omegaL*exp_sigma*(xL*dsigma_dy - yL*dsigma_dx); // line vanishes in axisymmetry
                   kxyL =  0.5*invlapse * exp_sigma * (yL*domega_dy - xL*domega_dx)
                          +0.5*invlapse * SQR(invrbar)*xL*yL*rns_omegaL*exp_sigma*( xL*dsigma_dy - yL*dsigma_dx )
                          +    invlapse * SQR(invrbar)*xL*yL*e2alp*( yL*dalpha_dx - xL*dalpha_dy );
                   kyyL = -    invlapse * xL * exp_sigma * domega_dy
                          +    invlapse * SQR(yL*invrbar)*rns_omegaL*e2alp*(yL*dalpha_dx-xL*dalpha_dy)
                          +0.5*invlapse * SQR(xL*invrbar)*rns_omegaL*exp_sigma*(yL*dsigma_dx-xL*dsigma_dy);
                   kxzL =  0.5*invlapse * yL * exp_sigma * domega_dz;                   
                   kyzL = -0.5*invlapse * xL * exp_sigma * domega_dz;
                   kzzL =      invlapse * rns_omegaL * e2alp * ( yL*dalpha_dx - xL*dalpha_dy ); // Vanishes in axisymmetry

                   betaxL =  rns_omegaL*yL;
                   betayL = -rns_omegaL*xL;
                   betazL =  0.;

                }


                if ( set_shift > 0 ) {
                   betax[idx] = betaxL; 
                   betay[idx] = betayL;
                   betaz[idx] = betazL;
                }

		if ( set_lapse > 0 ) {
		     alp[idx] = alpL;
		}

		// Fill Hydro velocities, changing to Valencia quantities and units where eos_kappa =/= 1
		if ( rns_velocity_sqL < 0 )
			rns_velocity_sqL = 0.;
		CCTK_REAL magv = sqrt( rns_velocity_sqL );
		//CCTK_REAL magv = -( rns_omegaL - *ns_omega ); // From T.S.'s conversation with Stergioulas
		CCTK_REAL velxL, velyL;
		if ( rbarL > 0 ) {
			velxL = -yL*magv*invrbar/sqrt(exp_sigma);
			velyL =  xL*magv*invrbar/sqrt(exp_sigma);
		} else {
                        velxL=0.;
                        velyL=0.;
                }
		CCTK_REAL velzL = 0;
		CCTK_REAL w_lorentzL = 1./sqrt( 1 - ( gxxL*SQR(velxL) + gyyL*SQR(velyL) + 2.*gxyL*velxL*velyL )); // Do we need to recreate v^2 with metric?
	
		// Fill Hydro densities and pressure
		CCTK_REAL pressL, rhoL, epsL;
		if ( rns_pressureL < pmin ) {
			// Outside star, set atmosphere with polytrope
			rhoL = rho_atmosphere;
			pressL = pmin; 
			epsL = epsmin; 
			velxL = 0.;
			velyL = 0.;
			velzL = 0.;
                        w_lorentzL = 1.;
		} else {

			if ( strcmp(eos_type,"poly")==0 ) {
		        	pressL = rns_pressureL; // Positive pressure
				rhoL = ( rns_energyL + pressL ) * exp(-rns_enthalpyL);
				//pressL = eos_Kappa * pow( rhoL, eos_Gamma ); // Make sure it's consistent
				//epsL = pressL/( rhoL*(eos_Gamma-1.)); 
				epsL = fmax( ( rns_energyL - rhoL )/rhoL , pressL/( rhoL*(eos_Gamma-1.) ) );
			} else {
                                
				int n_nearest = n_eosentries/2;
				rhoL = n0_at_e(rns_energyL, log_n0_tab, log_e_tab, n_eosentries,
                	               &n_nearest)*MB*KSCALE*SQ(C);
				epsL = ( rns_energyL - rhoL )/rhoL;
                                pressL = rns_pressureL;
			}
                       if ( add_pert ) {

                          /* Get r_surf(muL) */
                          /* In future, create an r_surf(muL) array to interpolate on outside 3D loop */
                          CCTK_REAL r_surface=0;
                          if ( *ns_oblateness != 1. ) {
                             /* Interpolate rsurf(mu) for muL */
                             if ( muL > (MDIV-2)*dmu || muL < 2*dmu ) {
                                /* Near pole or equator hack */
                                int app_idx = (floor)(muL/dmu);
                                r_surface = rsurf[app_idx];
                             } else {
                                CCTK_INT n_nearest=(int)( muL/dmu );
                                r_surface = rsurf[n_nearest]; 
                             }
                          } else {
                            r_surface = *ns_radius_iso;
                          }
                          if ( isnan(rL/r_surface) || r_surface < 0.001 ) {
                             #pragma omp critical
	                     CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                                  "At r=%g, mu=%g, muL=%g, r_surface=%g. ns_radius_iso=%g.", rL, muL, muL2, r_surface, *ns_radius_iso);
                          }
                          CCTK_REAL clean_rho=rhoL;
                          perturb_star(pert_idx, rL/r_surface, muL, xL, yL, muL2, &rhoL, &epsL, &pressL );
                          #pragma omp critical
                          {
                          fprintf( out_pert, "%g %g %g %.19g %.19g %.19g %.19g %.19g %.19g %.19g\n", xL, yL, zL, muL, muL2, rL, r_surface, rL/r_surface, 
                                   clean_rho, rhoL-clean_rho );
                          }
                       }
		}



		// Redundant check for physicality and atmosphere again
		if ( rhoL < rho_atmosphere || pressL < pmin ) {
			rhoL = rho_atmosphere;
			pressL = pmin; 
			epsL = epsmin; 
			velxL = 0.;
			velyL = 0.;
			velzL = 0.;
                        w_lorentzL = 1.;
		}

		if ( (pressL <= 0 || epsL < 0 || w_lorentzL < 1 || velxL > 1 || velyL > 1 || isnan(pressL) 
		     || isnan(epsL) || isnan(w_lorentzL) || isnan(velxL) || isnan(velyL) || isnan(alpL)
		     || isnan(betaxL) || isnan(betayL) || alpL <= 0 )) {
                        #pragma omp critical
			CCTK_INFO("Will be aborting the setting of initial data due to unphysical primitives or metric.");
			CCTK_VInfo(CCTK_THORNSTRING,"     Location:     (x,y,z) = (%g,%g,%g), r=%g, (s,mu) = (%g,%g) with r_e=%g, e(r_e)=%g.",xL,yL,zL,rL,sL,muL,*ns_radius_iso,e_surface);
			CCTK_VInfo(CCTK_THORNSTRING,"     Primitives:   rho=%g, P=%g, eps=%g, W=%g, v_x=%g, v_y=%g",rhoL,pressL,epsL,w_lorentzL,velxL,velyL);
			CCTK_VInfo(CCTK_THORNSTRING,"     Spacetime:    g=(%g,%g,%g,%g,%g,%g), alp=%g, beta=(%g,%g,0).",
				gxxL,gxyL,gxzL,gyyL,gyzL,gzzL, alpL, betaxL, betayL);
			CCTK_VInfo(CCTK_THORNSTRING,"     Interpolated Metric Potentials: alpha=%g, gamma=%g, rho=%g, omega=%g",rns_alphaL, rns_gamaL, rns_rhoL, rns_omegaL);
			CCTK_VInfo(CCTK_THORNSTRING,"     Interpolated Hydro: energy=%g, P=%g, enthalpy=%g, v^2=%g",rns_energyL, rns_pressureL, exp(rns_enthalpyL), rns_velocity_sqL);
                        CCTK_INT idx_center=IDX2D( nghosts+1, nghosts+1 );
			CCTK_VInfo(CCTK_THORNSTRING,"     NS Center Check: energy=%g, pressure=%g, enthalpy=%g, v^2=%g",rnsid_energy[idx_center],
							rnsid_pressure[idx_center],rnsid_enthalpy[idx_center],rnsid_velocity_sq[idx_center]);
			CCTK_WARN(0,"Error in setting RNS initial data.");

		}

		// Fill GFs
		gxx[idx] = gxxL;
                gxy[idx] = gxyL;
                gxz[idx] = gxzL;
		gyy[idx] = gyyL;
		gyz[idx] = gyzL;
		gzz[idx] = gzzL;

		kxx[idx] = kxxL;
		kxy[idx] = kxyL;
		kxz[idx] = kxzL;
		kyy[idx] = kyyL;
		kyz[idx] = kyzL;
		kzz[idx] = kzzL;

		rho[idx] = rhoL;
		eps[idx] = epsL;
		press[idx] = pressL;
#ifdef HAVE_MAYA_WHISKY
		velx[idx] = velxL;
		vely[idx] = velyL;
		velz[idx] = velzL;
#else
                const CCTK_INT xidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 0);
                const CCTK_INT yidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 1);
                const CCTK_INT zidx = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 2);
                vel[xidx] = velxL;
                vel[yidx] = velyL;
                vel[zidx] = velzL;
                w_lorentz[idx] = w_lorentzL;
#endif

	
  } 
  CCTK_ENDLOOP3_ALL(setRNSID);
  if ( add_pert ) 
     fclose( out_pert );

  /* Free space */
  free(rns_alpha);
  free(rns_omega);
  free(rns_rho);
  free(rns_gama);

  free(rns_e);
  free(rns_press);
  free(rns_h);
  free(rns_v2);

  free(dalpha_dmu);
  free(domega_dmu);
  free(drho_dmu);
  free(dgama_dmu);

  free(dalpha_ds);
  free(domega_ds);
  free(drho_ds);
  free(dgama_ds);

}

/*******************/
/* Local functions */
/*******************/
void RNS_output_state( double r_ratio, double central_energy_density, double mass, double mass0, double Re,
		double Omega, double Omega_K, double J, int with_header)
{

   double I_45;
   if (Omega == 0) I_45=0.0;
   else I_45 = J/(Omega*1.0e45);

   if ( with_header > 0 ) 
	CCTK_INFO("Current NS parameters: Oblateness, energy_central, mass(Msun)(grav), Mass(Msun)(rest), Omega, Omega_K, I_45, J");

   CCTK_VInfo(CCTK_THORNSTRING,"Current NS parameters: %4.3f \t%4.1f \t%4.3f \t%4.3f \t%4.3f \t%4.1f \t%5.1f \t%4.2f \t%4.3f",
	r_ratio, central_energy_density, mass/MSUN, mass0/MSUN, Re/1.0e5, Omega, Omega_K, I_45, (C*J/(G*SQR(mass))) );

}

void RNS_output_polystate( double r_ratio, double central_energy_density, double mass, double mass0, double rns_radius, 
  double Omega, double Omega_K, double J, double eos_Gamma, double eos_Kappa, int with_header)
{
   double I;
   if ( Omega == 0 ) I = 0.0;
   else I = J/Omega;

  // Solver solves in K=1 units. Output in K =/= 1 units
  // Base relevant conversion factor is K^n
  double eos_n = 1./(eos_Gamma-1.);
  double kapn = pow( eos_Kappa, eos_n );
  double ikapn = 1./kapn;

   if ( with_header > 0 ) 
	CCTK_INFO("Current NS parameters: Oblateness, energy_central, mass(grav), Mass(rest), rns_radius, Omega, Omega_K, I_45, J");

   CCTK_VInfo(CCTK_THORNSTRING,"Current NS parameters: %4.3g \t%4.3g \t%4.3g \t%4.3g \t%4.3g \t%4.3g \t%5.3g \t%4.2g \t%4.3g",
	r_ratio, ikapn*central_energy_density, sqrt(kapn)*mass, sqrt(kapn)*mass0, sqrt(kapn)*rns_radius, sqrt(ikapn)*Omega, sqrt(ikapn)*Omega_K, I, J/(SQR(mass0)) );

}

static CCTK_INT RNS_input_solution( CCTK_REAL *rho, CCTK_REAL *gama, CCTK_REAL *alpha, CCTK_REAL *omega, CCTK_REAL *energy, 
                                    CCTK_REAL *pressure, CCTK_REAL *enthalpy, CCTK_REAL *velocity_sq, 
                                    CCTK_REAL *ns_radius, CCTK_REAL *ns_radius_iso, CCTK_REAL *ns_mass, CCTK_REAL *ns_oblateness, CCTK_REAL *ns_omega )
{
  DECLARE_CCTK_PARAMETERS;

  char *filename;
  char format_str[2048];
  char line[2048];
  FILE *file;

  int npts_intrad = radial_spokes+2*nghosts;
  int npts_intang = angular_spokes+2*nghosts;
  int max_fname_chars = strlen(out_dir) + 1 + strlen(rns_file_in) + 1;

  filename =  (char *) malloc ( max_fname_chars ); // Maximum string size.
  assert(filename);

  /* Does file have a path structure? If so, assume full path has been given. */
  int has_pathstructure = 0;
  char check_filename[max_fname_chars+1];
  sprintf( check_filename, filename );
  char * fname_char = strchr(rns_file_in,'/');
  if ( fname_char != NULL )
        has_pathstructure=1;

  if ( has_pathstructure ) {
     sprintf(filename,"%s",rns_file_in);
  } else {
     sprintf(filename,"%s/%s", out_dir, rns_file_in);
  }

  if ( verbose )
	CCTK_VInfo(CCTK_THORNSTRING,"Reading in a previous RNS solution from file %s using format string %%lf.",filename);

  file = fopen(filename,"r");
  if (!file) {
	CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"RNS_input_solution: Could not open file %s for reading.",filename);
	return -1;
  }

  CCTK_INT nghosts_in = ( rns_file_ghostzone_size <0 ? nghosts : rns_file_ghostzone_size);
  CCTK_INT idx_off = ( rns_file_index_starts_at_one ? 1 : 0 );
  sprintf(format_str,"%%d\t%%d\t%%lf\t%%lf\t%%lf\t%%lf\t%%lf\t%%lf\t%%lf\t%%lf\t%%lf\t%%lf");

  CCTK_REAL eos_n = 1./(eos_Gamma-1.);
  CCTK_REAL kapn = pow( eos_Kappa, eos_n );
  CCTK_REAL ikapn = 1./pow( eos_Kappa, eos_n );

  int nlines=0, dlines=0;
  while ( fgets( line, sizeof(line), file ) != NULL ) {

	nlines++;
        if ( line[0] == '#' && nlines==2 ) {
           sscanf( line, "# ns_radius = %lf (Iso), %lf (Sch). ns_mass = %lf, ns_oblateness = %lf, ns_omega = %lf", 
                   ns_radius_iso, ns_radius, ns_mass, ns_oblateness, ns_omega);
           if ( rns_file_has_polytropic_units ) {
              *ns_radius_iso *= sqrt(kapn);
              *ns_radius *= sqrt(kapn);
              *ns_mass *= sqrt(kapn);
           }

        } 
	if ( line[0] == '#' ) continue; /* Ignore other header lines */

	if ( !strcmp( line, "\n" )) continue; /* We've come across a blank line */

	int iL,jL, idxL;
	CCTK_REAL sL, muL, muL2, rhoL, gamaL, alpL, omegaL, energyL, pressL, enthL, velsqL;

	/* Data line */
	sscanf( line, format_str, &iL, &jL, &sL, &muL, &muL2, &rhoL, &gamaL, &alpL, &omegaL, &energyL, &pressL, &enthL, &velsqL);

        if ( rns_file_has_polytropic_units ) {
           omegaL *= sqrt(ikapn);
           pressL *= ikapn;
           energyL *= ikapn;
        }
        
        idxL = IDX2D( (iL-idx_off) + (nghosts-nghosts_in), (jL-idx_off) + (nghosts-nghosts_in) );
	rho[idxL] = rhoL;
	gama[idxL] = gamaL;
	alpha[idxL] = alpL;
	omega[idxL] = omegaL;
	energy[idxL] = energyL;
	pressure[idxL] = pressL;
	enthalpy[idxL] = enthL;
	velocity_sq[idxL] = velsqL;
	dlines++;

  }
  if ( ! feof(file) ) {
     CCTK_INFO("Did not reach end of file.");
     strerror(errno);
  }

  /* Fill in ghostzones! */ 
  if ( nghosts != nghosts_in ) {
     fill_GA_ghostzones( rho, gama, alpha, omega, 
        energy, pressure, enthalpy, velocity_sq );
  }

  int lines_expected = (radial_spokes+2*nghosts_in)*(angular_spokes+2*nghosts_in);
  if ( dlines != lines_expected ) {
	double set_nghosts = 0.25*( -(radial_spokes + angular_spokes) + sqrt( SQR(radial_spokes) + SQR(angular_spokes) - 2*radial_spokes*angular_spokes + 4*dlines ) );
	CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Read %d data lines (of %d lines), which does not equal %d.  Please set nghosts to %g.",dlines,nlines,lines_expected,set_nghosts);
        fclose(file);
        free(filename);
	return -1;
  } else {
        fclose(file);
        free(filename);
        return 1;
  }

}

static CCTK_INT RNS_output_solution( CCTK_REAL *rho, CCTK_REAL *gama, CCTK_REAL *alpha, CCTK_REAL *omega, CCTK_REAL *energy, 
   CCTK_REAL *pressure, CCTK_REAL *enthalpy, CCTK_REAL *velocity_sq, 
   CCTK_REAL ns_radius, CCTK_REAL ns_radius_iso, CCTK_REAL ns_mass, CCTK_REAL ns_oblateness, CCTK_REAL ns_omega )
{
  DECLARE_CCTK_PARAMETERS;

  char *filename;
  char format_str_real[2048]; // XXX fixed size
  FILE *file;

  filename =  (char *) malloc ( strlen(out_dir) + 1 + strlen(rns_file) + 1 );
  assert(filename);
  sprintf(filename,"%s/%s", out_dir, rns_file);
  if ( verbose )
	CCTK_VInfo(CCTK_THORNSTRING,"Writing out solution of RNS to %s",filename);

  file = fopen(filename,"w");
  if (!file) {
	CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"RNS_output_solution: Could not open file %s",filename);
	return -1;
  }

  /* Header */
  fprintf(file,"# RNSID Solution: eps_central=%g, eos_info=%s (N=%g if applicable)\n",
	  central_energy_density, eos_info, 1./(eos_Gamma-1));
  fprintf(file,"# ns_radius = %lf (Iso), %lf (Sch). ns_mass = %lf, ns_oblateness = %lf, ns_omega = %lf\n", ns_radius_iso, 
          ns_radius, ns_mass, ns_oblateness, ns_omega);
  fprintf(file,"# 1:i 2:j 3:s 4:mu\n");
  fprintf(file,"#    Metric Potentials:  5:rho 6:gamma 7:alpha 8:omega\n");
  fprintf(file,"#    Hydro Solution:     9:energy(rho0+rho0*eps) 10:pressure 11:enthalpy 12:v^2\n");

  sprintf(format_str_real,"%%d\t%%d\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n",
           out_format,out_format,out_format,out_format,out_format,out_format,out_format,out_format,
	   out_format,out_format);

  CCTK_REAL ds = SMAX / (SDIV - 1.0);
  CCTK_REAL dmu = 1. / (MDIV - 1.0);
  /* Includes ghost zones */
  int npts_intrad = radial_spokes+2*nghosts;
  int npts_intang = angular_spokes+2*nghosts;
  for ( int i=0; i<npts_intrad; i++ ) {
     for ( int j=0; j<npts_intang; j++ ) {

	CCTK_REAL sL = (i-nghosts)*ds;
	CCTK_REAL muL = (j-nghosts)*dmu;
        CCTK_REAL muL2 = (j-nghosts)*dmu;
        CCTK_INT idx = IDX2D(i,j); 

	fprintf(file,format_str_real,i,j,sL,muL,muL2,
		rho[idx],gama[idx],alpha[idx],omega[idx],
		energy[idx],pressure[idx],enthalpy[idx],velocity_sq[idx]);

	// Add empty line between j's
	if ( j==npts_intang-1 )
		fprintf(file,"\n");

    }
  }
  fclose(file);
  free(filename);

  return 1;

}

static CCTK_INT perturb_star( const CCTK_INT pert_idx, const CCTK_REAL r_rs, const CCTK_REAL mu, const CCTK_REAL xx, const CCTK_REAL yy, CCTK_REAL r_inv, CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *press )
{

  DECLARE_CCTK_PARAMETERS;

  // Get central density from central e, assuming ideal gas and consistent pressure
  CCTK_REAL central_rho = central_energy_density - *press/( eos_Gamma - 1. );

  // Density perturbation applied. Need density from pressure and EoS
  // Assume polytrope here.

  CCTK_REAL drho=0.;
  CCTK_REAL dpress=0.;
  if ( pert_idx == 1 ) {
    // drho = pert_amp * central_rho * cos( 0.5*PI*r_rs );
     //dpress = eos_Gamma * (*press)* ( drho/ (*rho) );
     drho = pert_amp * central_rho*(SQR(r_rs))*(SQR(xx)-SQR(yy))*SQR(r_inv);
	 
  }
  if ( pert_idx == 2 ) {
      
	     drho = pert_amp * central_rho * sin( PI*r_rs )*(SQR(xx)-SQR(yy))*SQR(r_inv);
//       drho = pert_amp * central_rho * SQR(r_rs)*(SQR(xx)-SQR(yy))*SQR(r_inv);       
	          
  }

  *rho = *rho + drho;
  // Make consistent
  *press = eos_Kappa * pow(*rho,eos_Gamma);
  *eps = eos_Kappa * pow(*rho,eos_Gamma-1) / (eos_Gamma - 1. ); 

  return 0;
}

void RNS_InitGA(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT nrad = radial_spokes+2*nghosts;
  CCTK_INT nang = angular_spokes+2*nghosts;

  /* Initial scalars */
  *ns_radius = -1;
  *ns_radius_iso = -1;
  *ns_mass = -1;
  *ns_oblateness = -1;

  /* Initialize Solution */
  #pragma omp parallel for
  for ( int idx=0; idx<(nrad*nang); idx++ ) {
         rnsid_rho[idx] = INITVAL;
         rnsid_gama[idx] = INITVAL; 
         rnsid_alpha[idx] = INITVAL;
         rnsid_omega[idx] = INITVAL;
         rnsid_energy[idx] = INITVAL; 
         rnsid_pressure[idx] = INITVAL;
         rnsid_enthalpy[idx] = INITVAL;
         rnsid_velocity_sq[idx] = INITVAL;
  }

}

static CCTK_REAL *rns_allocate_array(CCTK_INT npoints, const char *name)
{
  DECLARE_CCTK_PARAMETERS;

  if (npoints<=0) {
    CCTK_WARN(0,"can't allocate array with npoints <=0");
  }
  if (name==NULL) {
    CCTK_WARN(0,"give a name");
  }

  CCTK_REAL *result;
  result=(CCTK_REAL *) malloc(sizeof(CCTK_REAL)*npoints);
  if (result==NULL) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "%s allocation for npoints=%d failed",name,npoints);
  }
  for (int i=0;i<npoints;i++) result[i]=0;
  const CCTK_REAL mbyt=npoints*sizeof(CCTK_REAL)/(1024.0*1024.0);
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"allocated array %s with %d elements -> %g MB",name,npoints,mbyt);
  }
  return result;
}

static void fill_GA_ghostzones( 
        CCTK_REAL * rnsid_rho, CCTK_REAL * rnsid_gama, CCTK_REAL * rnsid_alpha, CCTK_REAL * rnsid_omega, 
        CCTK_REAL * rnsid_energy, CCTK_REAL * rnsid_pressure, CCTK_REAL * rnsid_enthalpy, 
        CCTK_REAL * rnsid_velocity_sq )
{

     DECLARE_CCTK_PARAMETERS;
     int npts_intrad = radial_spokes+2*nghosts;

     /* Fill the sides */
     int muidx_pole = angular_spokes+nghosts-1; // index of highest physical mu
     int sidx_inf = radial_spokes+nghosts-1;  // index of highest physical s 
     for ( int j=0; j<=nghosts; j++ )   /* Angular ghostzones */
         for ( int i=nghosts-1; i<=sidx_inf; i++ ) { 
             CCTK_INT idx_gz, idx_eq, idx_pole;
            /* Inner ghosts are across equator: purely symmetric */
           idx_gz=IDX2D(i,j);
           idx_eq=IDX2D(i,2*nghosts-j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_eq];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_eq];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_eq];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_eq];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_eq];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_eq];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_eq];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_eq];

            /* Outer ghosts are across the pole */
            idx_gz=IDX2D(i, muidx_pole+j);
            idx_pole=IDX2D(i, muidx_pole-j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_pole];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_pole];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_pole];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_pole];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_pole];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_pole];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_pole];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_pole];

     }
     for ( int i=0; i<=nghosts; i++ ) /* Radial ghostzones */
         for ( int j=nghosts-1; j<=muidx_pole; j++ ){
             CCTK_INT idx_gz, idx_center, idx_inf;
            /* Inner ghosts in radial direction are purely symmetric */
            idx_gz=IDX2D(i,j);
            idx_center=IDX2D( 2*nghosts-i, j );
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_center];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_center];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_center];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_center];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_center];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_center];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_center];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_center];

            /* Outer ghosts in radial direction set to constant of last value (infinity is infinity) */
            idx_gz=IDX2D( sidx_inf+i, j );
            idx_inf=IDX2D( sidx_inf, j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_inf];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_inf];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_inf];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_inf];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_inf];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_inf];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_inf];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_inf];

     }

     /* Fill the corners */
     for ( int j=0; j<=nghosts; j++ ){
         for ( int i=0; i<=nghosts; i++ ) {

             CCTK_INT idx_gz, idx_c1, idx_c2, idx_c3, idx_c4;

            /* Low-mu, low-s corner purely symmetry */
            idx_gz=IDX2D(i,j);
            idx_c1=IDX2D(2*nghosts-i,2*nghosts-j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_c1];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_c1];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_c1];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_c1];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_c1];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_c1];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_c1];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_c1];
      
            /* Low-mu, high-s corner */
            idx_gz=IDX2D(sidx_inf+i,j);
            idx_c2=IDX2D(sidx_inf,2*nghosts-j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_c2];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_c2];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_c2];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_c2];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_c2];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_c2];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_c2];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_c2];

            /* Low-s, high-mu corner */
            idx_gz=IDX2D(i,muidx_pole+j);
            idx_c3=IDX2D(2*nghosts-i,muidx_pole-j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_c3];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_c3];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_c3];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_c3];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_c3];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_c3];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_c3];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_c3];

            /* High-s, high-mu corner */
            idx_gz=IDX2D(sidx_inf+i,muidx_pole+j);
            idx_c4=IDX2D(sidx_inf,muidx_pole-j);
            rnsid_rho[idx_gz]                  = rnsid_rho[idx_c4];
            rnsid_gama[idx_gz]                 = rnsid_gama[idx_c4];
            rnsid_alpha[idx_gz]                = rnsid_alpha[idx_c4];
            rnsid_omega[idx_gz]                = rnsid_omega[idx_c4];
            rnsid_energy[idx_gz]               = rnsid_energy[idx_c4];
            rnsid_pressure[idx_gz]             = rnsid_pressure[idx_c4];
            rnsid_enthalpy[idx_gz]             = rnsid_enthalpy[idx_c4];
            rnsid_velocity_sq[idx_gz]          = rnsid_velocity_sq[idx_c4];

         }
     }

}

