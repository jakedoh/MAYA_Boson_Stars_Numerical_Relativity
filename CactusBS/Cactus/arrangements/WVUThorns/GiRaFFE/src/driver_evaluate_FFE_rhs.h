#ifndef DRIVER_EVALUATE_MHD_RHS_H_
#define DRIVER_EVALUATE_MHD_RHS_H_

/* PRIVATE FUNCTIONS, Called within driver_evaluate_MHD_rhs.C ONLY */
static void reconstruct_set_of_prims_PPM_GRFFE(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
                                               const gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l, CCTK_REAL *temporary);

static void compute_TUPmunu
(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,const CCTK_REAL *dX,CCTK_REAL **metric,const gf_and_gz_struct *prims,
 const CCTK_REAL *gupxy,const CCTK_REAL *gupxz,const CCTK_REAL *gupyz,
 CCTK_REAL **TUPmunu);

static void A_i_rhs_no_gauge_terms(const int A_dirn,const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,
                                   CCTK_REAL *phi_interped,CCTK_REAL *cmax_1,CCTK_REAL *cmin_1,CCTK_REAL *cmax_2,CCTK_REAL *cmin_2, CCTK_REAL *A3_rhs);

static void Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,const CCTK_REAL *dX,CCTK_REAL **in_vars,const CCTK_REAL *psi6phi,
                                                           CCTK_REAL *shiftx_iphjphkph,CCTK_REAL *shifty_iphjphkph,CCTK_REAL *shiftz_iphjphkph,
                                                           CCTK_REAL *alpha_iphjphkph,CCTK_REAL *alpha_Phi_minus_betaj_A_j_iphjphkph,CCTK_REAL *alpha_sqrtg_Ax_interp,
                                                           CCTK_REAL *alpha_sqrtg_Ay_interp,CCTK_REAL *alpha_sqrtg_Az_interp,
                                                           CCTK_REAL *psi6phi_rhs,CCTK_REAL *Ax_rhs,CCTK_REAL *Ay_rhs,CCTK_REAL *Az_rhs);

static void add_fluxes_and_source_terms_to_hydro_rhss(const int flux_dirn,const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,const CCTK_REAL *dX,
                                                      CCTK_REAL **metric,const gf_and_gz_struct *in_prims,CCTK_REAL **TUPmunu,
                                                      const int numvars_reconstructed,const gf_and_gz_struct *out_prims_r,const gf_and_gz_struct *out_prims_l,
                                                      CCTK_REAL *cmax,CCTK_REAL *cmin,
                                                      CCTK_REAL *st_x_flux,CCTK_REAL *st_y_flux,CCTK_REAL *st_z_flux,
                                                      CCTK_REAL *st_x_rhs,CCTK_REAL *st_y_rhs,CCTK_REAL *st_z_rhs);

#endif /* DRIVER_EVALUATE_MHD_RHS_H_ */
