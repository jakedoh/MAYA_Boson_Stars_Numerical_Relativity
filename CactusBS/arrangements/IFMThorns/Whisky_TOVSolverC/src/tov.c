/* file    tov.c
 * author  Frank Loeffler, converted from fortran thorn by Ian Hawke
 * date    2002/10/21
 * desc    TOV initial data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "EinsteinBase/Constants/src/constants.h"

#include "tov.h"

CCTK_REAL * TOV_Surface=0;
CCTK_REAL * TOV_R_Surface=0;
CCTK_REAL * TOV_Atmosphere=0;

CCTK_REAL * TOV_r_1d=0;
CCTK_REAL * TOV_rbar_1d=0;
CCTK_REAL * TOV_press_1d=0;
CCTK_REAL * TOV_phi_1d=0;
CCTK_REAL * TOV_m_1d=0;

#include "utils.inc"

void TOV_C_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (TOV_Solve_for_TOVs != 3)
    if (TOV_Solve_for_TOVs == 2)
    {
      CCTK_WARN(1, "TOV_Solve_for_TOVs is depreciated. "
                   "Use TOV_Enforce_Interpolation=\"yes\" instead.\n");
      if (CCTK_ParameterSet("TOV_Enforce_Interpolation",
                            "Whisky_TOVSolverC",
                            "true"))
        CCTK_WARN(0, "Error while steering this parameter.\n");
      else
        CCTK_WARN(1, "Steered this parameter for now accordingly.\n");
    }
    else
      CCTK_WARN(1, "TOV_Solve_for_TOVs is depreciated. "
                   "Use TOV_Enforce_Interpolation instead.\n");
}

/* centered differencing with one-sided differencing at the boundary */
#define DIFF_X(a) (((i==0)?(a[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i  , j, k)]): \
                    (i==(cctk_lsh[0]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i  , j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]))/\
                   CCTK_DELTA_SPACE(0))
#define DIFF_Y(a) (((j==0)?(a[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j  , k)]): \
                    (j==(cctk_lsh[1]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i, j  , k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]))/\
                    CCTK_DELTA_SPACE(1))
#define DIFF_Z(a) (((k==0)?(a[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k  )]): \
                    (k==(cctk_lsh[2]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i, j, k  )] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]))/\
                    CCTK_DELTA_SPACE(2))


/*@@
   @routine    TOV_Source_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler - converted fortran routine by Ian Hawke
   @desc
      The source terms for the ODEs. These are equations (2), (3), (4)
      and (18) from the Baumgarte notes.
      That is the vector in order is (P, m, phi, rbar).
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV_C_Source_RHS(CCTK_REAL r, CCTK_REAL K, CCTK_REAL Gamma,
                      CCTK_REAL old_data[4], CCTK_REAL source_data[4])
{
  CCTK_REAL LOCAL_TINY, PI;
  CCTK_REAL press, rho, eps, mu, m, phi;
  CCTK_REAL log_rbar_over_r, r_minus_two_m;

  LOCAL_TINY = 1.0e-35;
  PI=4.0*atan(1.0);

  press           = old_data[0];
  if (press < LOCAL_TINY)
    press = LOCAL_TINY;
  m               = old_data[1];
  phi             = old_data[2];
  log_rbar_over_r = old_data[3];

  rho = pow(press / K, 1.0 / Gamma);
  eps = press / (Gamma - 1.0) / rho;
  mu  = rho * (1.0 + eps);

  r_minus_two_m = r - 2.0 * m;

  if ((r==0.0) && (m==0.0))
  {
    /* source_data[0] = 0.0; */
    source_data[1] = 0.0;
    source_data[2] = 0.0;
    source_data[3] = 0.0;
  }
  else
  {
    source_data[2] = (m + 4*PI * r*r*r * press) / r_minus_two_m / r;
    /* source_data[0] = -(press + mu) * source_data[2]; */
    source_data[0] = -(press + mu) *
                     (m + 4*PI * r*r*r * press) / r_minus_two_m / r;
    source_data[1] = 4*PI * r*r * mu;
    source_data[3] = (sqrt(r) - sqrt(r_minus_two_m)) / r / sqrt(r_minus_two_m);
  }
}

/*@@
   @routine    TOV_Integrate_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
      Integrates the ODEs using RK4.
      We rescale at the end to match to a Schwarzschild exterior.
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/

void TOV_C_Integrate_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL LOCAL_TINY;

  CCTK_INT star, star_i, i, TOV_Surface_Index;
  CCTK_REAL old_data[4], source_data[4], in_data[4], new_data[4],
            k1[4], k2[4], k3[4], k4[4];
  CCTK_REAL Surface_Mass, factor, local_rho;
  CCTK_REAL rho_central;

  LOCAL_TINY = 1.0e-20;

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);
  assert(TOV_Atmosphere!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  /* do it for all stars */
  for (star=0; star < TOV_Num_TOVs; star++)
  {
    /* remember array index */
    star_i = star * TOV_Num_Radial;

    /* check for parameters */
    if ((TOV_Rho_Central[star]==0.0) && (whisky_rho_central>0.0))
      rho_central=whisky_rho_central;
    else
      rho_central=TOV_Rho_Central[star];

    if (rho_abs_min > 0.0)
      TOV_Atmosphere[star] = rho_abs_min;
    else
      TOV_Atmosphere[star] = rho_central * rho_rel_min;

    if (initial_rho_abs_min > 0.0)
      TOV_Atmosphere[star] = initial_rho_abs_min;
    else
      if (initial_rho_rel_min > 0.0)
        TOV_Atmosphere[star] = rho_central * initial_rho_rel_min;

    if (initial_atmosphere_factor > 0.0)
      TOV_Atmosphere[star] *= initial_atmosphere_factor;

    /* Set conformal state like set in parameter file if we do not use
     * the old initial data. In this case we have to use what we get */
    if (!TOV_Use_Old_Initial_Data)

    /* clear arrays first */
    TOV_C_fill(&(TOV_press_1d[star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_m_1d    [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_phi_1d  [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_rbar_1d [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_r_1d    [star_i]), TOV_Num_Radial, 0.0);

    /* set start values */
    TOV_press_1d[star_i] = TOV_K[star] *
                            pow(rho_central, TOV_Gamma[star]);
    /* TOV_r_1d    [star_i] = LOCAL_TINY;
    TOV_rbar_1d [star_i] = LOCAL_TINY;*/

    /* build TOV_r_1d[] */
    for (i=star_i+1; i < star_i+TOV_Num_Radial; i++)
      TOV_r_1d[i] = TOV_r_1d[i-1] + TOV_dr[star];

    TOV_Surface[star] = -1.0;

    /* loop over all radii */
    for (i=star_i; (i < star_i+TOV_Num_Radial-1) &&
                            (TOV_Surface[star] < 0.0); i++)
    {
      old_data[0] = TOV_press_1d[i];
      old_data[1] = TOV_m_1d[i];
      old_data[2] = TOV_phi_1d[i];
      if (TOV_rbar_1d[i]==TOV_r_1d[i])
        old_data[3] = 0.0;
      else
        old_data[3] = log(TOV_rbar_1d[i] / TOV_r_1d[i]);

      in_data[0] = old_data[0];
      in_data[1] = old_data[1];
      in_data[2] = old_data[2];
      in_data[3] = old_data[3];

      TOV_C_fill(source_data, 4, 0.0);

      TOV_C_Source_RHS(TOV_r_1d[i],
                     TOV_K[star], TOV_Gamma[star],
                     in_data, source_data);
      k1[0] = TOV_dr[star] * source_data[0];
      k1[1] = TOV_dr[star] * source_data[1];
      k1[2] = TOV_dr[star] * source_data[2];
      k1[3] = TOV_dr[star] * source_data[3];
      in_data[0] = old_data[0] + 0.5 * k1[0];
      in_data[1] = old_data[1] + 0.5 * k1[1];
      in_data[2] = old_data[2] + 0.5 * k1[2];
      in_data[3] = old_data[3] + 0.5 * k1[3];

      TOV_C_Source_RHS(TOV_r_1d[i]+ 0.5 * TOV_dr[star],
                       TOV_K[star], TOV_Gamma[star],
                       in_data, source_data);
      k2[0] = TOV_dr[star] * source_data[0];
      k2[1] = TOV_dr[star] * source_data[1];
      k2[2] = TOV_dr[star] * source_data[2];
      k2[3] = TOV_dr[star] * source_data[3];
      in_data[0] = old_data[0] + 0.5 * k2[0];
      in_data[1] = old_data[1] + 0.5 * k2[1];
      in_data[2] = old_data[2] + 0.5 * k2[2];
      in_data[3] = old_data[3] + 0.5 * k2[3];

      TOV_C_Source_RHS(TOV_r_1d[i]+ 0.5 * TOV_dr[star],
                       TOV_K[star], TOV_Gamma[star],
                       in_data, source_data);
      k3[0] = TOV_dr[star] * source_data[0];
      k3[1] = TOV_dr[star] * source_data[1];
      k3[2] = TOV_dr[star] * source_data[2];
      k3[3] = TOV_dr[star] * source_data[3];
      in_data[0] = old_data[0] + k3[0];
      in_data[1] = old_data[1] + k3[1];
      in_data[2] = old_data[2] + k3[2];
      in_data[3] = old_data[3] + k3[3];

      TOV_C_Source_RHS(TOV_r_1d[i]+ TOV_dr[star],
                       TOV_K[star], TOV_Gamma[star],
                       in_data, source_data);
      k4[0] = TOV_dr[star] * source_data[0];
      k4[1] = TOV_dr[star] * source_data[1];
      k4[2] = TOV_dr[star] * source_data[2];
      k4[3] = TOV_dr[star] * source_data[3];

      new_data[0] = old_data[0] + (k1[0] + k4[0] + 2.0 * (k2[0] + k3[0])) /6.0;
      new_data[1] = old_data[1] + (k1[1] + k4[1] + 2.0 * (k2[1] + k3[1])) /6.0;
      new_data[2] = old_data[2] + (k1[2] + k4[2] + 2.0 * (k2[2] + k3[2])) /6.0;
      new_data[3] = old_data[3] + (k1[3] + k4[3] + 2.0 * (k2[3] + k3[3])) /6.0;

      TOV_press_1d[i+1] = new_data[0];
      TOV_m_1d    [i+1] = new_data[1];
      TOV_phi_1d  [i+1] = new_data[2];
      TOV_rbar_1d [i+1] = TOV_r_1d[i+1] * exp(new_data[3]);

      /* otherwise the code crashes later */
      if (TOV_press_1d[i+1] < 0.0)
          TOV_press_1d[i+1] = 0.0;

      local_rho = pow(TOV_press_1d[i+1] / TOV_K[star], 1.0 / TOV_Gamma[star]);

      /* scan for the surface */
      /*  if ( (local_rho <= TOV_Atmosphere[star]) ||*/
      if ( (local_rho <= 0.0) ||
           (TOV_press_1d[i+1] <= 0.0) )
      {
        TOV_Surface[star]   = TOV_r_1d[i];
        TOV_R_Surface[star] = TOV_rbar_1d[i];
        TOV_Surface_Index = i;
      }
    }
    if (TOV_Surface[star] < 0.0)
      CCTK_WARN(0, "Did not integrate out to surface of the star! "
                   "Increase TOV_dr or TOV_Num_Radial and rerun");

    Surface_Mass = TOV_m_1d[TOV_Surface_Index];
    factor = 0.5 * (sqrt(TOV_Surface[star] *
                         (TOV_Surface[star] - 2.00 * Surface_Mass)) +
                    TOV_Surface[star] - Surface_Mass) /
                    TOV_rbar_1d[TOV_Surface_Index];

    TOV_R_Surface[star] *= factor;
    for (i=star_i; i < star_i+TOV_Num_Radial; i++)
    {
      TOV_rbar_1d[i] *= factor;
      TOV_phi_1d[i]  -= TOV_phi_1d[TOV_Surface_Index] -
                        0.5 * log(1.0 - 2.0 * Surface_Mass / TOV_Surface[star]);
      /* match to Schwarzschield */
      if (i > TOV_Surface_Index)
      {
        TOV_press_1d[i] = TOV_K[star] * pow(TOV_Atmosphere[star],
                                            TOV_Gamma[star]);
        TOV_rbar_1d [i] = 0.5 *
                          (sqrt(TOV_r_1d[i]*(TOV_r_1d[i] - 2.0*Surface_Mass)) +
                          TOV_r_1d[i] - Surface_Mass);
        TOV_m_1d[i]     = Surface_Mass;
        TOV_phi_1d[i]   = 0.5 * log( 1.0 - 2.0 * Surface_Mass / TOV_r_1d[i]);
      }
    }
  }
  CCTK_INFO("Integrated TOV equation");
  /* do some info */
  CCTK_VInfo(CCTK_THORNSTRING, "Information about the TOVs used:");
  CCTK_VInfo("", "TOV    radius     mass  mass(g) cent.rho rho(cgi)        K   K(cgi)    Gamma");
  for (i=0; i<TOV_Num_TOVs; i++)
    if (TOV_Gamma[i]==2.0)
      CCTK_VInfo("","  %d  %8g %8g %8.3g %8g %8.3g %8g %8.3g %8g",
                 i+1, TOV_R_Surface[i],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1]*CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i]/pow(CONSTANT_G_cgi,3.0)/
                                    pow(CONSTANT_Msolar_cgi,2.0)*
                                    pow(CONSTANT_c_cgi,6.0),
                 TOV_K[i],
                 TOV_K[i]*pow(CONSTANT_G_cgi,3.0)*
                          pow(CONSTANT_Msolar_cgi,2.0)/
                          pow(CONSTANT_c_cgi,4.0),
                 TOV_Gamma[i]);
    else
      CCTK_VInfo("","  %d  %8g %8g %8.3g %8g %8.3g %8g %8g",
                 i+1, TOV_R_Surface[i],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1]*CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i]/pow(CONSTANT_G_cgi,3.0)/
                                    pow(CONSTANT_Msolar_cgi,2.0)*
                                    pow(CONSTANT_c_cgi,6.0),
                 TOV_K[i], TOV_Gamma[i]);

}

/*----------------------------------------------------------------------------*/

/* utility routine
 * recursive search-routine for arrays
 * here used to look for the last index in an ordered array with its
 * value < goal
 */
CCTK_INT TOV_C_find_index(CCTK_INT   array_size,
                          CCTK_REAL *array,
                          CCTK_REAL  goal,
                          CCTK_INT   lower_index,
                          CCTK_INT   upper_index)
{
  CCTK_INT middle_index;

  if (lower_index >= (upper_index-1))
      return lower_index;

  middle_index = (lower_index + upper_index) /2;

  if (array[middle_index] < goal)
    return TOV_C_find_index(array_size, array, goal, middle_index, upper_index);
  else
    return TOV_C_find_index(array_size, array, goal, lower_index, middle_index);
}

/* utility rountine
 * interpolates from (thorn-internal) 1D-data to Cactus 3D-grid */
/* input is all but *press_point *phi_point and *r_point */
void TOV_C_interp_tov_isotropic(
                                CCTK_INT  star,
                                CCTK_REAL *TOV_press_1d,
                                CCTK_REAL *TOV_phi_1d,
                                CCTK_REAL *TOV_rbar_1d,
                                CCTK_REAL *TOV_r_1d,
                                CCTK_REAL *r,
                                CCTK_REAL surface,
                                CCTK_REAL *press_point,
                                CCTK_REAL *phi_point,
                                CCTK_REAL *r_point,
                                CCTK_INT  TOV_warn_small_grid,
                                CCTK_INT  TOV_zero_atmosphere,
                                CCTK_REAL *TOV_Atmosphere)
{
  DECLARE_CCTK_PARAMETERS
  CCTK_INT  left_index;
  CCTK_REAL h, M;

  if (*r < 0.0)
    CCTK_WARN(0, "Negative radius found");
  if (*r < TOV_rbar_1d[1])
    *r=TOV_rbar_1d[1];
  if (*r > TOV_rbar_1d[TOV_Num_Radial-2])
    if (TOV_warn_small_grid)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Grid not large enough - last point: %.20f.",
                 TOV_rbar_1d[TOV_Num_Radial-1]);
    else
    {
      if (!TOV_zero_atmosphere)
        *press_point= TOV_K[star] * pow(TOV_Atmosphere[star],
                                        TOV_Gamma[star]);
      else
        *press_point=0.0;
      M = 0.5 * TOV_r_1d[TOV_Num_Radial-1] *
                (1.0 - exp(2.0*TOV_phi_1d[TOV_Num_Radial-1]));
      *r_point=(2* *r+M)*(2* *r+M)*0.25/ *r;
      *phi_point=0.5*log(1-2*M/ *r_point);
      return;
    }

  if (TOV_Fast_Interpolation)
    left_index = TOV_C_find_index(TOV_Num_Radial-1, TOV_rbar_1d, *r, 0,
                                  TOV_Num_Radial-1);
  else
  {
    left_index=0;
    while( (left_index < TOV_Num_Radial-2) &&
           (TOV_rbar_1d[left_index+1] < *r) )
      left_index++;
  }

  h = (*r - TOV_rbar_1d[left_index]) /
      (TOV_rbar_1d[left_index+1] - TOV_rbar_1d[left_index]);
  *r_point     = (1.0 - h) * TOV_r_1d[left_index] +
                        h  * TOV_r_1d[left_index+1];
  *phi_point   = (1.0 - h) * TOV_phi_1d[left_index] +
                        h  * TOV_phi_1d[left_index+1];
  if (!TOV_zero_atmosphere ||
      (*r_point < surface))
    *press_point = (1.0 - h) * TOV_press_1d[left_index] +
                          h  * TOV_press_1d[left_index+1];
  else
    *press_point = 0.0;
}

/*@@
   @routine    TOV_Exact
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
       Schedule routine for interpolation of 1D to 3D grid
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV_C_Exact(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL *press_point, *rho_point, *eps_point,
            *mu_point, *phi_point, *r_point;

  CCTK_INT  LSH_MAX_I;

  CCTK_INT i,j,k, i3D, star, star_i;
  CCTK_REAL *r_to_star;
  CCTK_REAL g_diag, max_g_diag, max_rho, det, sqrt_det;
  CCTK_REAL my_velx, my_vely, my_velz, my_psi4;
  CCTK_REAL vlowx, vlowy, vlowz;
  CCTK_REAL D_h_w, PI, local_tiny;

  CCTK_INT tov_lapse;

  tov_lapse = CCTK_EQUALS(initial_lapse, "tov");

  PI=4.0*atan(1.0);
  local_tiny=1.0e-14;
  /* remember index of last member of array */
  LSH_MAX_I = CCTK_GFINDEX3D(cctkGH,
                             cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1);

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);
  assert(TOV_Atmosphere!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  /* allocate local arrays */
  r_to_star   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  press_point = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  rho_point   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  eps_point   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  mu_point    = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  phi_point   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  r_point     = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));

  /* clear initial data */
  if (TOV_Clear_Initial_Data > 0 && !(TOV_Use_Old_Initial_Data))
  {
    TOV_C_fill(kxx, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kxy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kxz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kyy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kyz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kzz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gxx, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gyy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gzz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gxy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gxz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gyz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(alp, LSH_MAX_I+1, 1.0);

    if (*shift_state != 0)
    {
      TOV_C_fill(betax, LSH_MAX_I+1, 0.0);
      TOV_C_fill(betay, LSH_MAX_I+1, 0.0);
      TOV_C_fill(betaz, LSH_MAX_I+1, 0.0);
    }

    if (!TOV_Use_Old_Matter_Initial_Data)
    {
      TOV_C_fill(rho,        LSH_MAX_I+1, 0.0);
      TOV_C_fill(dens,       LSH_MAX_I+1, 0.0);
      TOV_C_fill(eps,        LSH_MAX_I+1, 0.0);
      TOV_C_fill(press,      LSH_MAX_I+1, 0.0);
      TOV_C_fill(tau,        LSH_MAX_I+1, 0.0);
      TOV_C_fill(w_lorentz,  LSH_MAX_I+1, 0.0);
      TOV_C_fill(sx,         LSH_MAX_I+1, 0.0);
      TOV_C_fill(sy,         LSH_MAX_I+1, 0.0);
      TOV_C_fill(sz,         LSH_MAX_I+1, 0.0);
      TOV_C_fill(velx,       LSH_MAX_I+1, 0.0);
      TOV_C_fill(vely,       LSH_MAX_I+1, 0.0);
      TOV_C_fill(velz,       LSH_MAX_I+1, 0.0);
    }
  }
  /* use the fast interpolation? only useful for testing this */
  if (TOV_Fast_Interpolation == 0)
      CCTK_INFO("Interpolating the slow way.");

  /* go over all 3D-grid points */
  for(i=0; i<cctk_lsh[0]; i++)
   for(j=0; j<cctk_lsh[1]; j++)
    for(k=0; k<cctk_lsh[2]; k++)
    {
      i3D=CCTK_GFINDEX3D(cctkGH, i, j, k);
      /* remember the old conformal factor to the power of 4 */
      my_psi4=1.0;

      for (star=0; star<TOV_Num_TOVs; star++)
      {
        r_to_star[star] =
          sqrt( (x[i3D]-TOV_Position_x[star]) *
                (x[i3D]-TOV_Position_x[star]) +
                (y[i3D]-TOV_Position_y[star]) *
                (y[i3D]-TOV_Position_y[star]) +
                (z[i3D]-TOV_Position_z[star]) *
                (z[i3D]-TOV_Position_z[star]) );
        star_i = star * TOV_Num_Radial;

        /* do the actual interpolation */
        TOV_C_interp_tov_isotropic(star,
                               &(TOV_press_1d[star_i]), &(TOV_phi_1d[star_i]),
                               &(TOV_rbar_1d[star_i]), &(TOV_r_1d[star_i]),
                               &(r_to_star[star]), TOV_Surface[star],
                               &(press_point[star]),
                               &(phi_point[star]), &(r_point[star]),
                               0, 0, TOV_Atmosphere);

        /* is some perturbation wanted? */
        if (Perturb[star] == 0)
          rho_point[star] = pow(press_point[star]/TOV_K[star],
                                1.0/TOV_Gamma[star]);
        else
          rho_point[star] = pow(press_point[star]/TOV_K[star],
                                1.0/TOV_Gamma[star]) *
                            (1.0 +
                             Pert_Amplitude[star] *
                               cos(PI/2.0 * r[i3D] / TOV_R_Surface[star]));

        eps_point[star] = press_point[star] / (TOV_Gamma[star] - 1.0)
                                            /  rho_point[star];
        mu_point[star]  = rho_point[star] * (1.0 + eps_point[star]);
      }
      /* find out from which star we want to have the data */
      if (CCTK_EQUALS(TOV_Combine_Method, "maximum"))
      {
        /* to do this, we use here simply the max of the gxx-value */
        star=0;
        max_g_diag = 0.0;
        max_rho = rho_point[0];
        for (star_i=0; star_i<TOV_Num_TOVs; star_i++)
        {
          g_diag = (r_point[star_i] / (r_to_star[star_i] + 1.0e-30)) *
                   (r_point[star_i] / (r_to_star[star_i] + 1.0e-30));
          if ((g_diag - max_g_diag) > local_tiny)
          {
            max_g_diag=g_diag;
            star=star_i;
          }
          if ((rho_point[star_i] - max_rho) > local_tiny)
          {
            max_rho=rho_point[star_i];
            star=star_i;
          }
        }
        /* handle initial data */
        if (TOV_Use_Old_Initial_Data)
        {
          /* check metric */
          if ((my_psi4 * gxx[i3D] < max_g_diag) &&
              (my_psi4 * gyy[i3D] < max_g_diag) &&
              (my_psi4 * gzz[i3D] < max_g_diag))
          {
            if (TOV_Conformal_Flat_Three_Metric)
            {
              psi[i3D] = pow(max_g_diag, 0.25);
              my_psi4 = max_g_diag;
            }
            else
            {
              gxx[i3D] = max_g_diag/my_psi4;
              gyy[i3D] = max_g_diag/my_psi4;
              gzz[i3D] = max_g_diag/my_psi4;
            }
          }
          /* check matter */
          if (TOV_Use_Old_Matter_Initial_Data)
          {
            if (rho[i3D] > max_rho)
            {
              /* we do not need this array element anymore, since we use
               * the initial data, so lets use it */
              star=0;
              max_rho          =rho[i3D];
              eps_point[star]  =eps[i3D];
              press_point[star]=press[i3D];
              my_velx=velx[i3D];
              my_vely=vely[i3D];
              my_velz=velz[i3D];
            }
            else
            {
              if (tov_lapse)
                alp[i3D] = exp(phi_point[star]);
              my_velx=TOV_Velocity_x[star];
              my_vely=TOV_Velocity_y[star];
              my_velz=TOV_Velocity_z[star];
            }
          }
          else
          {
            if (tov_lapse)
              alp[i3D] = exp(phi_point[star]);
            my_velx=TOV_Velocity_x[star];
            my_vely=TOV_Velocity_y[star];
            my_velz=TOV_Velocity_z[star];
          }
        }
        else /* do not use old initial data */
        {
          /* no psi, since it is 1.0 here */
          /* but maybe we want to have it != 1.0 */
          if (TOV_Conformal_Flat_Three_Metric)
          {
            psi[i3D] = pow(max_g_diag, 0.25);
            my_psi4 = max_g_diag;
            gxx[i3D] = gyy[i3D] = gzz[i3D] = 1.0;
          }
          else
          {
            gxx[i3D] = max_g_diag;
            gyy[i3D] = max_g_diag;
            gzz[i3D] = max_g_diag;
          }
          if (tov_lapse)
            alp[i3D] = exp(phi_point[star]);
          my_velx=TOV_Velocity_x[star];
          my_vely=TOV_Velocity_y[star];
          my_velz=TOV_Velocity_z[star];
        }

        /* set to defined velocity. default is 0.0 because other velocities
         * violate Einsteins equations */
        velx[i3D] = my_velx;
        vely[i3D] = my_vely;
        velz[i3D] = my_velz;

        w_lorentz[i3D] = 1/sqrt(1.0-(
                                  gxx[i3D] * velx[i3D] * velx[i3D]+
                                  gyy[i3D] * vely[i3D] * vely[i3D]+
                                  gzz[i3D] * velz[i3D] * velz[i3D]+
                                2*gxy[i3D] * velx[i3D] * vely[i3D]+
                                2*gxz[i3D] * velx[i3D] * velz[i3D]+
                                2*gyz[i3D] * vely[i3D] * velz[i3D])*
                                my_psi4);

        rho[i3D] = max_rho;
        eps[i3D] = eps_point[star];
        press[i3D] = press_point[star];

      }
      else if (CCTK_EQUALS(TOV_Combine_Method, "average"))
      {
        /* here we 'average' all values in a more intelligent way */
        if (TOV_Use_Old_Matter_Initial_Data)
          max_rho=rho[i3D];
        else
        {
          max_rho=0.0;
          rho[i3D] = 0.0;
        }
        star=-1;
        for (star_i=0; star_i<TOV_Num_TOVs; star_i++)
        {
          if (tov_lapse)
            alp[i3D] *= exp(phi_point[star_i]);
          if (TOV_Conformal_Flat_Three_Metric)
          {
            /* This is a hack, since it does not check, if the input data is
             * really conformally flat. It simply assumes this by only using
             * gxx */
            my_psi4 = (r_point[star_i] * r_point[star_i] /
                       (r_to_star[star_i] * r_to_star[star_i] + 1.0e-30)) /
                      my_psi4 + pow(psi[i3D], 4.0) * gxx[i3D];
            psi[i3D] = pow(my_psi4, 0.25);
            if (!TOV_Use_Old_Initial_Data)
              gxx[i3D] = gyy[i3D] = gzz[i3D] = 1.0;
          }
          else
            gxx[i3D] += (r_point[star_i] * r_point[star_i] /
                         (r_to_star[star_i] * r_to_star[star_i] + 1.0e-30)) /
                        my_psi4;
          rho[i3D] += rho_point[star_i];
          eps[i3D] += eps_point[star_i];
          press[i3D] += press_point[star_i];
          /* we still have to know if we are inside one star - and which */
          if ( (rho_point[star_i] > max_rho) &&
               (rho_point[star_i] > TOV_Atmosphere[star_i]))
          {
            max_rho=rho_point[star_i];
            star=star_i;
          }

        if(rho[i3D] <= TOV_Atmosphere[star_i])
            rho[i3D] = TOV_Atmosphere[star_i];

        }

        if (TOV_Conformal_Flat_Three_Metric)
        {
          my_psi4 -= ((TOV_Num_TOVs+TOV_Use_Old_Initial_Data-1)/my_psi4);
          psi[i3D] = pow(my_psi4, 0.25);
        }
        else
        {
          gxx[i3D] -= ((TOV_Num_TOVs+TOV_Use_Old_Initial_Data-1)/my_psi4);
          gyy[i3D] = gxx[i3D];
          gzz[i3D] = gxx[i3D];
        }

        /* set to defined velocity. default is 0.0 because other velocities
         * violate the constraints */
        if (star > -1)
        {
          velx[i3D] = TOV_Velocity_x[star];
          vely[i3D] = TOV_Velocity_y[star];
          velz[i3D] = TOV_Velocity_z[star];
        }

        w_lorentz[i3D] = 1/sqrt(1.0-(
                                  gxx[i3D] * velx[i3D] * velx[i3D]+
                                  gyy[i3D] * vely[i3D] * vely[i3D]+
                                  gzz[i3D] * velz[i3D] * velz[i3D]+
                                2*gxy[i3D] * velx[i3D] * vely[i3D]+
                                2*gxz[i3D] * velx[i3D] * velz[i3D]+
                                2*gyz[i3D] * vely[i3D] * velz[i3D]) * my_psi4);
      }
      /* now all combine methods are doing the same */
      if ((TOV_Num_TOVs==1) && (rho[i3D] < TOV_Atmosphere[0]))
        velx[i3D] = vely[i3D] = velz[i3D] = 0.0;

      SpatialDet(gxx[i3D],gxy[i3D],gxz[i3D],gyy[i3D],gyz[i3D],gzz[i3D],&det);

      sqrt_det=sqrt(det) * pow(my_psi4, 1.5);
      dens[i3D] = sqrt_det * w_lorentz[i3D] * rho[i3D];

      /* this variable is only for temporal storage */
      D_h_w = dens[i3D] * (1 + eps[i3D] + press[i3D]/rho[i3D]) * w_lorentz[i3D];

      vlowx = gxx[i3D]*velx[i3D] + gxy[i3D]*vely[i3D] + gxz[i3D]*velz[i3D];
      vlowy = gxy[i3D]*velx[i3D] + gyy[i3D]*vely[i3D] + gyz[i3D]*velz[i3D];
      vlowz = gxz[i3D]*velx[i3D] + gyz[i3D]*vely[i3D] + gzz[i3D]*velz[i3D];

      sx[i3D] = D_h_w * vlowx;
      sy[i3D] = D_h_w * vlowy;
      sz[i3D] = D_h_w * vlowz;

      /* One could use D_h_w here, but it would introduce more error */
      tau[i3D] = sqrt_det * (rho[i3D]  *w_lorentz[i3D]*w_lorentz[i3D] *
                                         (1.0 + eps[i3D]) +
                              press[i3D]*(w_lorentz[i3D]*w_lorentz[i3D]-1.0) ) -
                 dens[i3D];
    }
  /* if used, recalculate the derivatives of the conformal factor */
  i3D = cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];
  switch(TOV_Populate_Timelevels)
  {
    case 3:
        TOV_Copy(i3D, gxx_p_p,  gxx);
        TOV_Copy(i3D, gyy_p_p,  gyy);
        TOV_Copy(i3D, gzz_p_p,  gzz);
        TOV_Copy(i3D, rho_p_p,  rho);
        TOV_Copy(i3D, eps_p_p,  eps);
        TOV_Copy(i3D, velx_p_p, velx);
        TOV_Copy(i3D, vely_p_p, vely);
        TOV_Copy(i3D, velz_p_p, velz);
        TOV_Copy(i3D, dens_p_p, dens);
        TOV_Copy(i3D, tau_p_p,  tau);
        TOV_Copy(i3D, sx_p_p,   sx);
        TOV_Copy(i3D, sy_p_p,   sy);
        TOV_Copy(i3D, sz_p_p,   sz);
        TOV_Copy(i3D, w_lorentz_p_p, w_lorentz);
    case 2:
        TOV_Copy(i3D, gxx_p,  gxx);
        TOV_Copy(i3D, gyy_p,  gyy);
        TOV_Copy(i3D, gzz_p,  gzz);
        TOV_Copy(i3D, rho_p,  rho);
        TOV_Copy(i3D, eps_p,  eps);
        TOV_Copy(i3D, velx_p, velx);
        TOV_Copy(i3D, vely_p, vely);
        TOV_Copy(i3D, velz_p, velz);
        TOV_Copy(i3D, dens_p, dens);
        TOV_Copy(i3D, tau_p,  tau);
        TOV_Copy(i3D, sx_p,   sx);
        TOV_Copy(i3D, sy_p,   sy);
        TOV_Copy(i3D, sz_p,   sz);
        TOV_Copy(i3D, w_lorentz_p, w_lorentz);
  }
  CCTK_INFO("Done interpolation.");

  /* free local arrays */
  free(r_to_star);
  free(press_point);
  free(rho_point);
  free(eps_point);
  free(mu_point);
  free(phi_point);
  free(r_point);
}

void TOV_Prepare_Fake_Evolution(CCTK_ARGUMENTS)
{
    if (CCTK_ParameterSet("TOV_Populate_Timelevels",
                          "Whisky_TOVSolverC",
                          "1"))
        CCTK_WARN(0,
                  "Could not prepare for fake evolution - steering failed\n");
}

#include "external.inc"
