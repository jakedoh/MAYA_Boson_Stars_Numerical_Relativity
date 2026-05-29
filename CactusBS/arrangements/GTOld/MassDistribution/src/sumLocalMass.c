#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define KRANC_C
#include "GenericFD.h"

#include "MassDistribution.h"

#define SQR(x) ((x) * (x))

/*************************/
/*****static helpers******/
/*************************/
static void get_interior_domain_specs(CCTK_ARGUMENTS, CCTK_INT imin[3], CCTK_INT imax[3]);
static void parse_per_volume_parameters(CCTK_ARGUMENTS, 
        CCTK_REAL *VolWeights[], 
        enum bin_mode my_bin_mode[], 
        enum bin_variable_type bin_vartype[], CCTK_REAL *bin_variable_ptr[],
        enum mass_variable_type mass_vartype[], CCTK_REAL *mass_variable_ptr[], 
        CCTK_REAL x0L[], CCTK_REAL y0L[], CCTK_REAL z0L[]);

/* if you feel like doing a good debugging run of this, then add the following
 * line to TOVSolver's tov.c somewhere after the ODE has been solved:
 *
    double mtot = 0.;
    #define SQR(x) ((x)*(x))
    #define CUBE(x) ((x)*(x)*(x))
    star = 0;
    for (i=0; i < TOV_Surface_Index+1; i++)
    {
        double r = TOV_r_1d[i];
        double rbar = TOV_rbar_1d[i];
        double P = TOV_press_1d[i];
        double rho = pow(P/TOV_K[star],1/TOV_Gamma[star]);
        double eps = P/((TOV_Gamma[star]-1)*rho);
        double M = TOV_m_1d[i];
        double mu = rho * (1 + eps);
        double dP = -(mu+P)*(M+4*M_PI*CUBE(r)*P)/(r*(r-2*M));
        double deps = 1/TOV_Gamma[star]*dP/rho;
        double drho = 1/TOV_Gamma[star]*rho/P*dP;

        double drbar = (1.+(sqrt(r)-sqrt(r-2*M))/sqrt(r-2*M))*rbar/r;
        double W = 1/sqrt(1-SQR(TOV_Velocity_z[star])*SQR(r/rbar));
        double dW = SQR(TOV_Velocity_z[star])*CUBE(W)*r/SQR(rbar)*(1-r/rbar*drbar);
        double E = (1+eps)*W+(SQR(W)-1)/W*P/rho;
        double dE = W*deps + (1+eps)*dW + (SQR(W)-1)/W*(dP/rho-P/SQR(rho)*drho) + P/rho*(dW+dW/SQR(W));
        double dM = 4*M_PI*rho*W * SQR(r)/sqrt(1-2*M/r); // rest mass only
        double dM_dE = dM/dE;
        if(r>0.)
        mtot += 0.5*4*M_PI*( 
                pow(TOV_press_1d[i]/TOV_K[star],1/TOV_Gamma[star])*SQR(TOV_r_1d[i])/sqrt(1-2*TOV_m_1d[i]/TOV_r_1d[i]) + 
                pow(TOV_press_1d[i+1]/TOV_K[star],1/TOV_Gamma[star])*SQR(TOV_r_1d[i+1])/sqrt(1-2*TOV_m_1d[i+1]/TOV_r_1d[i+1])
                ) * W * TOV_dr[star];

        if(CCTK_MyProc(NULL)==0)
          fprintf(stderr, "%g %g %g %g\n", E, dM_dE, dE, r);
    }
    fprintf(stderr, "Total mass: %f, TOV_Mass: %f\n", mtot, TOV_m_1d[TOV_Surface_Index]);
 *
 * then run the code on the test parameter file and do in gnuplot:
 * plot [:] [0:] "matterlibs_tovtest/stderr" u ($1):(-$2), "matterlibs_tovtest/MassDistribution0.asc" u ($4):($5) w lp
 * the two lines should lie on top of each other (up to noise).
*/

/* from GenericFD: obtain range of indices to cover the interior domain (plus
 * the outer physical boundary) */
static void get_interior_domain_specs(CCTK_ARGUMENTS, CCTK_INT imin[3], CCTK_INT imax[3])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Summation limits */
  CCTK_INT   is_symbnd[6], is_physbnd[6], is_ipbnd[6];
  CCTK_INT   nboundaryzones[6], is_internal[6], is_staggered[6], shiftout[6];
  CCTK_INT   npoints;

  GenericFD_GetBoundaryInfo(cctkGH, cctk_ash, cctk_lsh, cctk_bbox,
                            cctk_nghostzones,
                            imin, imax, is_symbnd, is_physbnd, is_ipbnd);

  if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
    int ierr = GetBoundarySpecification(6, nboundaryzones, is_internal, is_staggered, shiftout);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
  } else {
    /* assume that nboundaryzones = nghostzones
     * leave is_internal, is_staggered, shiftout undefined
     */
    for (int dir = 0; dir<3; dir++)
    {
      for (int face = 0; face < 2; face++)
      {
         CCTK_INT indx = dir*2 + face;
         nboundaryzones[indx] = cctk_nghostzones[dir];
      }
    }
  }

  for (int dir = 0; dir<3; dir++)
  {
    for (int face = 0; face < 2; face++)
    {

       CCTK_INT indx = dir*2 + face;
       if (is_physbnd[indx]) {
   	   npoints=0;
       } else if (is_symbnd[indx]) {
	   npoints=nboundaryzones[indx];
       } else {
  	   npoints=cctk_nghostzones[dir];
       }

       switch(face)
       {
       case 0: /*lower*/
          imin[dir]=npoints;
	  break;
       case 1: /*upper*/
          imax[dir]=cctk_lsh[dir] - npoints;
	  break;
       default:
          CCTK_WARN(0,"internal error");
          return; /* NOTREACHED */
       }

    }
  }

  if (verbose >= 3) { 
	CCTK_VInfo(CCTK_THORNSTRING,"Extent: (%d,%d,%d) to (%d,%d,%d)",imin[0],imin[1],imin[2],imax[0],imax[1],imax[2]);
  }
}

static void parse_per_volume_parameters(CCTK_ARGUMENTS, 
        CCTK_REAL *VolWeights[], 
        enum bin_mode my_bin_mode[], 
        enum bin_variable_type bin_vartype[], CCTK_REAL *bin_variable_ptr[],
        enum mass_variable_type mass_vartype[], CCTK_REAL *mass_variable_ptr[], 
        CCTK_REAL x0L[], CCTK_REAL y0L[], CCTK_REAL z0L[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vartype;

  for(int vol = 0 ; vol < intvolumes ; vol++)
  {
    /* Get Weight Pointers */
    VolWeights[vol] = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_variable[vol]));
    assert(VolWeights[vol]);

    /* Get type (and pointer) to binning variable */
    MassDistribution_resolve_special_variable(cctkGH, bin_variable[vol], 
        &vartype, &bin_variable_ptr[vol], 
        NUM_SPECIAL_VARIABLE_TYPES, special_variables);
    bin_vartype[vol] = (enum bin_variable_type)vartype;
    assert(bin_vartype[vol] != binvar_gridfunction || 
           (bin_vartype[vol] == binvar_gridfunction && bin_variable_ptr[vol] != NULL));

    /* Get type (and pointer) to mass variable */
    MassDistribution_resolve_special_variable(cctkGH, mass_variable[vol], 
        &vartype, &mass_variable_ptr[vol], 
        NUM_MASS_VARIABLE_TYPES, mass_variables);
    mass_vartype[vol] = (enum mass_variable_type)vartype;
    assert(mass_vartype[vol] != massvar_gridfunction || 
           (mass_vartype[vol] == massvar_gridfunction && mass_variable_ptr[vol] != NULL));

    /* is the bining range explicitely given or should I find minima and maxima? */
    if(CCTK_Equals(bin_mode[vol], "automatic")) {
      my_bin_mode[vol] = automatic;
    } else if(CCTK_Equals(bin_mode[vol], "automatic-volume")) {
      my_bin_mode[vol] = automatic_volume;
    } else if(CCTK_Equals(bin_mode[vol], "manual")) {
      my_bin_mode[vol] = manual;
    } else if(CCTK_Equals(bin_mode[vol], "mass-fraction")) {
      my_bin_mode[vol] = by_mass_fraction;
    } else {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Unknown range mode '%s' for volume %d.", bin_mode[vol], vol);
      return; /* NOTREACHED */
    }
    assert(my_bin_mode[vol] <= by_mass_fraction);

    /* set up local origins for individual integration volumes */
    MassDistribution_get_centre(CCTK_PASS_CTOC, vol, &x0L[vol], &y0L[vol], &z0L[vol]);
  }
}

/*************************/
/****external interface***/
/*************************/

void MassDistribution_sumLocalMass(CCTK_ARGUMENTS, enum operation opcode)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rhoL,epsL,pressL,vxL,vyL,vzL,w2L, wL;
  CCTK_REAL gxxL,gyyL,gzzL,gxyL,gxzL,gyzL,detg;
  CCTK_REAL vlowxL,vlowyL,vlowzL;
  CCTK_REAL alphaL, betaxL, betayL, betazL;
  CCTK_REAL weightL,MassL;
  CCTK_REAL dx,dy,dz,dV;
  CCTK_REAL xL,yL,zL,r2L;
  CCTK_REAL x0L[intvolumes],y0L[intvolumes],z0L[intvolumes];
  CCTK_INT check_location[intvolumes];
  CCTK_REAL *VolWeights[intvolumes], *bin_variable_ptr[intvolumes];
  CCTK_REAL bin_variable_value;
  CCTK_REAL mass_variable_value, *mass_variable_ptr[intvolumes];
  CCTK_INT imin[3], imax[3];
  double symsigns[8][3];
  int ncopies;
  enum bin_mode my_bin_mode[intvolumes];
  enum bin_variable_type bin_vartype[intvolumes];
  enum mass_variable_type mass_vartype[intvolumes];
  int vol, bin, bin_upper, bin_lower;
  int num_negative_mass_points;

  /* check input */
  if(opcode > integrate)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, 
            "Unknown operation %d to sumLocalRestMass", (int)opcode);
    assert(0); /* NOTREACHED */
  }
  
  /* get symmetry description */
  query_symmetries(applied_symmetries, &ncopies, symsigns);

  /* Parse per volume parameters to get pointers and enums from strings */
  parse_per_volume_parameters(CCTK_PASS_CTOC, VolWeights, my_bin_mode,
      bin_vartype, bin_variable_ptr, mass_vartype, mass_variable_ptr, 
      x0L, y0L, z0L);
  for(int v = 0 ; v < intvolumes ; v++)
  {
    if(opcode == integrate)
      check_location[v] = 1;
    else if(opcode == find_minmax && my_bin_mode[v] == automatic_volume)
      check_location[v] = 1;
    else if(opcode == find_minmax && my_bin_mode[v] == by_mass_fraction)
      check_location[v] = 1;
    else
      check_location[v] = 0;
  }

  /* Loop over everything except ghost zones */ 
  get_interior_domain_specs(CCTK_PASS_CTOC, imin, imax);
  num_negative_mass_points = 0;
  for(int k = imin[2] ; k < imax[2] ; k++)
  {
    for(int j = imin[1] ; j < imax[1] ; j++)
    {
      for(int i = imin[0] ; i < imax[0] ; i++)
      {
        for(int s = 0 ; s < ncopies ; s++)
        {

          /* Information from the grid */
          int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
          
          /* anything that is too low is treated as vacuum and does not contribute to the mass */
          if ( rho[idx] <= rho_min )
            continue;
          
          rhoL = rho[idx];
          epsL = eps[idx];
          pressL = press[idx];
          
          xL = symsigns[s][0]*x[idx];
          yL = symsigns[s][1]*y[idx];
          zL = symsigns[s][2]*z[idx];
          
          gxxL = gxx[idx]; /* signs always cancel */
          gxyL = symsigns[s][0]*symsigns[s][1]*gxy[idx];
          gxzL = symsigns[s][0]*symsigns[s][2]*gxz[idx];
          gyyL = gyy[idx]; /* signs always cancel */
          gyzL = symsigns[s][1]*symsigns[s][2]*gyz[idx];
          gzzL = gzz[idx]; /* signs always cancel */
          
          vxL = symsigns[s][0]*velx[idx];
          vyL = symsigns[s][1]*vely[idx];
          vzL = symsigns[s][2]*velz[idx];

          alphaL = alp[idx];
          betaxL = symsigns[s][0]*betax[idx];
          betayL = symsigns[s][1]*betay[idx];
          betazL = symsigns[s][2]*betaz[idx];

          vlowxL = gxxL*vxL + gxyL*vyL + gxzL*vzL;
          vlowyL = gxyL*vxL + gyyL*vyL + gyzL*vzL;
          vlowzL = gxzL*vxL + gyzL*vyL + gzzL*vzL;

          w2L = 1. / ( 1. - (vlowxL*vxL + vlowyL*vyL + vlowzL*vzL) );
          wL = sqrt(w2L);

          /* Get volume element */
          dx = CCTK_DELTA_SPACE(0);
          dy = CCTK_DELTA_SPACE(1);
          dz = CCTK_DELTA_SPACE(2);
          dV = dx*dy*dz;
          
          detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
          
          for (vol=0; vol<intvolumes; vol++)
          {
            if ( !HistogramActive[vol] || !HistogramValid[vol] )
              continue;
          
            weightL = VolWeights[vol][idx];
            if ( weightL == 0 ) 
              continue;
          
            r2L = SQR(xL-x0L[vol]) + SQR(yL-y0L[vol]) + SQR(zL-z0L[vol]);
            if ( check_location[vol] && (r2L > SQR(int_wi_rad[vol])) )
              continue;
          
            /* get variable value */
            switch(bin_vartype[vol])
            {
              double v2SR, w2SR; /* temp vars for binvar_ESR */
              case binvar_gridfunction:
                bin_variable_value = bin_variable_ptr[vol][idx];
                break;
              case binvar_specific_energy:
                bin_variable_value = (1+epsL)*wL + (w2L - 1.)/wL*pressL/rhoL;
                break;
              case binvar_u0:
                bin_variable_value = wL*(alphaL - (betaxL*vlowxL+betayL*vlowyL+betazL*vlowzL));
                break;
              case binvar_ENewton:
                bin_variable_value = 0.5*(SQR(vxL) + SQR(vyL) + SQR(vzL)) - central_mass[vol]/sqrt(r2L); 
                break;
              case binvar_ESR:
                v2SR = SQR(vxL) + SQR(vyL) + SQR(vzL);
                w2SR = 1./(1.-v2SR);
                bin_variable_value =  0.5*w2SR*v2SR - central_mass[vol]/sqrt(r2L)*(1-SQR(w2SR-1));
                break;
              default:
                assert(0 && "internal error");
            }

            if ( bin_variable_minimum[vol] != bin_variable_maximum[vol] &&
                 (bin_variable_value < bin_variable_minimum[vol] ||
                  bin_variable_value >=  bin_variable_maximum[vol]) )
              continue;
                
            /* finally do something */
            switch(opcode)
            {
              case integrate:
                /* calculate the amount of rest mass in this cell */
                switch(mass_vartype[vol])
                {
                  case massvar_gridfunction:
                    mass_variable_value = mass_variable_ptr[vol][idx];
                    break;
                  case massvar_relativistic_rest_mass_density:
                    mass_variable_value = wL * rhoL;
                    break;
                  case massvar_debug1:
                    mass_variable_value = 0.123 * xL * xL;
                    break;
                  default:
                    assert(0 && "internal error");
                }
                MassL = sqrt(detg) * mass_variable_value * dV;
                if (MassL < 0)
                {
                  num_negative_mass_points += 1;
                  if (verbose >= 4) 
                    CCTK_VInfo(CCTK_THORNSTRING,"Negative MassL found in MassDistribution at (%g,%g,%g).  MassL=%g.", xL, yL, zL, MassL);
                }


                /* check for validity and set to ignore data */
                if(ignore_invalid_cells && 
                    (!isfinite(bin_variable_value) || !isfinite(MassL)))
                  continue;

                /* sum if within range of valid values */
                if(bin_variable_value >= HistogramMinimum[vol] && bin_variable_value < HistogramMaximum[vol])
                {
                  /* find bin number (bisection) */
                  for (bin_lower = 0, bin_upper = number_of_bins, bin = (bin_lower+bin_upper)/2 ; 
                       bin != bin_lower ;
                       bin = (bin_lower+bin_upper)/2)
                  {
                    if (bin_variable_value >= HistogramLabels[vol*NUMBER_OF_LABELS + bin])
                      bin_lower = bin;
                    else
                      bin_upper = bin;
                  }
                  assert(bin >= 0 && bin < number_of_bins);
                  MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + bin] += MassL*weightL;
                }

                /* accumulate to total mass */
                MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + TOTAL_MASS_BIN] += MassL*weightL;
                /* accumulate mass above and below threshold */
                if(bin_variable_value <= bin_variable_threshold[vol])
                  MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + BELOW_MASS_BIN] += MassL*weightL;
                else
                  MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + ABOVE_MASS_BIN] += MassL*weightL;

                break;
              case find_minmax:
                /* check for validity and set to ignore data */
                if(ignore_invalid_cells && !isfinite(bin_variable_value))
                  continue;

                /* check if this is bigger/smaller than the old max/min respectively */
                if(bin_variable_value < HistogramMinimum[vol])
                  HistogramMinimum[vol] = bin_variable_value;
                if(bin_variable_value > HistogramMaximum[vol])
                  HistogramMaximum[vol] = bin_variable_value;

                break;
              default:
                assert(0 && "internal error");
                break;
            } /* switch(opcode) */
          } /* for (vol=0; vol<intvolumes; vol++) */
        } /* for s */
      } /* for i */
    } /* for j */
  } /* for k */

  if ( num_negative_mass_points>0 ) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Local mass integrand in MassDistribution is negative in %d points.  This implies there's something wrong.", 
               num_negative_mass_points);
  }

  if(verbose >= 3)
  {
    for(vol = 0 ; vol < intvolumes ; vol++)
    {
      if(HistogramActive[vol] && HistogramValid[vol]) 
      {
        switch(opcode)
        {
          case integrate:
            CCTK_VInfo(CCTK_THORNSTRING,"Partial mass in volume %d is %f.", vol, 
                MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + TOTAL_MASS_BIN]);
            break;
          case find_minmax:
            if (my_bin_mode[vol] == automatic || my_bin_mode[vol] == automatic_volume || 
                my_bin_mode[vol] == by_mass_fraction)
            {
              CCTK_VInfo(CCTK_THORNSTRING,"Range found in volume %d is [%f,%f].", 
                      vol, HistogramMinimum[vol], HistogramMaximum[vol]);
            }
            break;
          default:
            assert(0 && "internal error");
        } /* switch opcode */
      } /* if HistogramActive */
    } /* for vol */
  } /* if verbose */

}
