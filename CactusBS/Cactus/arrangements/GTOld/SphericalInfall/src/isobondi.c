#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "util_String.h"

#define INITVALUE (42)
#define M_PI     3.14159265358979323846  /* pi */
#define PSIMAX   1e4

/* Scheduled routines */
void SphericalInfall_IsoBondi_Solve(CCTK_ARGUMENTS);
void SphericalInfall_IsoBondi_Populate(CCTK_ARGUMENTS);
void SphericalInfall_IsoBondi_WriteSolns(CCTK_ARGUMENTS);

/* Static functions */
static CCTK_INT calc_BondiSolution( CCTK_REAL r, CCTK_REAL * rho, CCTK_REAL * u, 
   const CCTK_REAL rSonic, const CCTK_REAL rhoSonic, 
   const CCTK_REAL M, const CCTK_REAL Mdot, const CCTK_REAL eos_K, 
   const CCTK_REAL Qdot );

// Inline helpers
static inline CCTK_REAL SQR (const CCTK_REAL x)
{
  return ((x)*(x));
}
static inline CCTK_REAL QUAD (const CCTK_REAL x)
{
  const CCTK_REAL y = ((x)*(x));
  return ((y)*(y));
}

void SphericalInfall_IsoBondi_Solve ( CCTK_ARGUMENTS ) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Initialize GridArray */
  for(int nsol=0; nsol<bondi_num_solutions; nsol++ ) {
     for( int n=0; n<bondi_maxnpts; n++ ) {
        int idx = bondi_maxnpts*nsol + n;
        r_bondi[idx] = INITVALUE;
        logr_bondi[idx] = INITVALUE;
        rho_bondi[idx] = INITVALUE;
        u_bondi[idx] = INITVALUE;
        int_energy_bondi[idx] = INITVALUE;
     }
  } 

  /* Solve all solutions */
  for( int nsol=0; nsol<bondi_num_solutions; nsol++ ) {

     /**********************/
     /* Process parameters */
     /**********************/
     CCTK_REAL M, M2, Mdot;
     CCTK_REAL rSonic, rmin, rmax, npts;

     M = bondi_central_mass[nsol];
     Mdot = bondi_mdot[nsol];
     M2 = SQR(M);

     rSonic = bondi_sonicpt_radius[nsol];
     rmin = M * bondi_rmin[nsol]; 
     rmax = M * bondi_rmax[nsol];
     npts = bondi_npts[nsol]; 

     /* Sound speed squared at sonic point */
     CCTK_REAL aSonic2, aSonic, uSonic2, uSonic;
     aSonic2 = M / ( 2*rSonic - 3*M );
     if ( aSonic2 > ( eos_gamma - 1. ) ) {
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "Problem with Bondi solution %d. aSonic^2 > gam-1.", nsol);
     }
     aSonic = sqrt(aSonic2); 
     uSonic2 = M / (2*rSonic);
     uSonic = sqrt(uSonic2);
    
     /* EoS K from solution at sonic point */
     CCTK_REAL rhoSonic, eos_K, Qdot, hSonic;
     hSonic = 1./( 1. - aSonic2/(eos_gamma-1.));
     rhoSonic = Mdot / ( 4. * M_PI * uSonic * SQR(rSonic) );
     eos_K = hSonic * aSonic2 * pow( rhoSonic, 1.-eos_gamma ) / eos_gamma;
     Qdot = SQR(hSonic) * ( 1. - 3.*uSonic2 );

     CCTK_VInfo(CCTK_THORNSTRING,"Bondi Solution %d:\t eos_K=%g, Mdot=%g, Qdot=%g",
                nsol, eos_K, Mdot, Qdot );
     CCTK_VInfo(CCTK_THORNSTRING,"      at Sonic Pt: r=%g, rho=%g, a=%g, h=%g, u=%g",
                rSonic, rhoSonic, aSonic, hSonic, uSonic );

     /* Set up minimum radius for solution */
     CCTK_REAL logrmin, dlogr;
     logrmin = log10(rmin);
     dlogr   = ( log10(rmax) - logrmin) / ( 1.*(npts-1) );

     /* Logarithmic spacing in r */
     CCTK_INT idx_Sonic = 0; 
     CCTK_REAL dist_rSonic = 1e30; 
     for ( int rn=0; rn < npts; rn++ ) {
         int idx = bondi_maxnpts*nsol + rn;
         logr_bondi[idx] = logrmin + dlogr*(rn-1);
         r_bondi[idx] = pow( 10, logr_bondi[idx] );
         if ( fabs(r_bondi[idx] - rSonic) < dist_rSonic ) {
            dist_rSonic = fabs(r_bondi[idx] - rSonic);
            idx_Sonic = idx;
         }
     }

     CCTK_REAL rhoL = rhoSonic;
     CCTK_REAL uL = uSonic;
     CCTK_INT  ierr;
     /* Start solving from sonic point outwards,                    */
     /*   We want the larger rho solution outside the sonic points */
     for ( int i=idx_Sonic; i<npts; i++ ) { 
         CCTK_REAL rSch = r_bondi[i];
         rhoL = 1.02*rhoL;
         ierr = calc_BondiSolution( rSch, &rhoL, &uL, rSonic, rhoSonic, M, Mdot, eos_K, Qdot );
         if ( rhoL < rho_min ) {
            rhoL = rho_min;
            uL   = eos_K * pow( rhoL, eos_gamma ) / ( eos_gamma - 1. );
         } 
         rho_bondi[i] = rhoL;
         u_bondi[i]   = uL;
         int_energy_bondi[i]   = eos_K * eos_gamma * pow( rhoL, eos_gamma ) / ( eos_gamma - 1.);
     }

     /* Solve from sonic point inwards                              */
     /*   We want the smaller rho solution outside the sonic points */
     rhoL = rhoSonic;
     uL = uSonic;
     for ( int i=idx_Sonic-1; i>=0; i-- ) {
         CCTK_REAL rSch = r_bondi[i];
         rhoL = 0.98*rhoL;
         ierr = calc_BondiSolution( rSch, &rhoL, &uL, rSonic, rhoSonic, M, Mdot, eos_K, Qdot );
         if ( rhoL < rho_min ) {
            rhoL = rho_min;
            uL   = eos_K * pow( rhoL, eos_gamma ) / ( eos_gamma - 1. );
         } 
         rho_bondi[i] = rhoL;
         u_bondi[i]   = uL;
         int_energy_bondi[i]   = eos_K * eos_gamma * pow( rhoL, eos_gamma ) / ( eos_gamma - 1.);
     }

  }

}

void SphericalInfall_IsoBondi_Populate( CCTK_ARGUMENTS ) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Populate with 0th solution */
#define NSOL 0
  const CCTK_REAL M = bondi_central_mass[NSOL]; 
  const CCTK_REAL rSonic = bondi_sonicpt_radius[NSOL];
  const CCTK_REAL Mdot = bondi_mdot[NSOL];
  const CCTK_REAL logrmin = log10(bondi_rmin[NSOL]);
  const CCTK_REAL dlogr   = ( log10(bondi_rmax[NSOL]) - logrmin) / ( 1.*(bondi_npts[NSOL]-1) );

  const CCTK_REAL aSonic2 = M / ( 2*rSonic - 3*M );
  const CCTK_REAL uSonic2 = M / (2*rSonic);
  const CCTK_REAL hSonic = 1./( 1. - aSonic2/(eos_gamma-1.));
  const CCTK_REAL rhoSonic = Mdot / ( 4. * M_PI * sqrt(uSonic2) * SQR(rSonic) );
  const CCTK_REAL eos_K = hSonic * aSonic2 * pow( rhoSonic, 1.-eos_gamma ) / eos_gamma;
  const CCTK_REAL Qdot = SQR(hSonic) * ( 1. - 3.*uSonic2 );

  /* Store the derivative of rho near rSch ~ 2.25 M */
  CCTK_REAL rminus, rplus, rho_m, rho_plus, drhodr, u_plus, u_minus;
  CCTK_INT jminus, jplus, ierr;
  
  rminus = 2.25 * M;
  jminus = floor( 0.5 + ( log10(rminus) - logrmin) / dlogr );
  rho_m = rho_bondi[jminus];
  u_minus   = u_bondi[jminus];
  ierr = calc_BondiSolution ( rminus, &rho_m, &u_minus, rSonic, rhoSonic,
         M, Mdot, eos_K, Qdot );
  
  rplus  = 0.25 * SQR(3.02) * M/1.01; // ~2.2575
  jplus = floor( 0.5 + ( log10(rplus) - logrmin) / dlogr );
  rho_plus = rho_bondi[jplus];
  u_plus   = u_bondi[jplus];
  ierr = calc_BondiSolution ( rplus, &rho_plus, &u_plus, rSonic, rhoSonic,
         M, Mdot, eos_K, Qdot );

  drhodr = 100*( rho_plus - rho_m )/M;

  const CCTK_INT set_lapse = CCTK_EQUALS(initial_lapse,"bondi_infall_iso");

  for( int k=0; k<cctk_lsh[2]; k++ ) {
     for( int j=0; j<cctk_lsh[1]; j++ ) {
        for(int i=0; i<cctk_lsh[0]; i++ ) {
 
           CCTK_INT idx=CCTK_GFINDEX3D(cctkGH, i, j, k);

           /* Assume coordinates are isotropic */
           CCTK_REAL xL, yL, zL, rIso;
           xL = x[idx];
           yL = y[idx];
           zL = z[idx];
           rIso = sqrt( SQR(xL) + SQR(yL) + SQR(zL) );

           /* Fill in geometry, smoothing over the puncture as in Faber et al. */
           /* Note: this is not constraint-satisfying                          */
           CCTK_REAL psi, psi4;
           if ( rIso < 0.5*M ) {
              //psi = 2.875 - 5.*SQR(rIso/M) + 6.*QUAD(rIso/M);
              psi = 2.875 - 5.*pow(rIso/M,2.) + 6.*pow(rIso/M,4.);
           } else {
              psi = fmin( 1.+ M/(2.*rIso), PSIMAX );
           }
           psi4 = QUAD(psi); 
           gxx[idx] = psi4;
           gyy[idx] = psi4;
           gzz[idx] = psi4;
           gxy[idx] = 0.;
           gxz[idx] = 0.;
           gyz[idx] = 0.;
           if ( set_lapse > 0 ) {
              alp[idx] = 1./SQR(psi);
           }

           kxx[idx] = 0.;
           kyy[idx] = 0.;
           kzz[idx] = 0.;
           kxy[idx] = 0.;
           kxz[idx] = 0.;
           kyz[idx] = 0.;

           /* Calculate the equivalent schwarzschild radius from the closest 
              solution point where r < r_i                                   */
           CCTK_REAL rSch, uSch, uIso, rhoL, epsL;
           if ( rIso > M ) {
              /* Far from BH */
              CCTK_INT rIdx;
              rSch = SQR( 2.*rIso + M )/( 4.*rIso);
              rIdx = floor( 0.5 + ( log10(rSch) - logrmin ) / dlogr );
              
              rhoL = rho_bondi[rIdx];
              uSch = u_bondi[rIdx];
              ierr = calc_BondiSolution ( rSch, &rhoL, &uSch, rSonic, rhoSonic,
                M, Mdot, eos_K, Qdot );

              uIso = uSch / ( 1. - SQR( M/(2.*rIso) ) );
           } else {
              if ( rIso > 0.5*M ) {
                 /* Use linearization of rho in vicinity of BH */
                 rhoL = rho_m + drhodr*rIso*(rIso-M)/M;
              } else {
                 /* rho within the horizon, to be excised later */
                 rhoL = ( rho_m-drhodr*M/4.)*0.5*( 1.-cos(2.*M_PI*rIso/M) );
              }
              uIso = u_minus*rIso/M;
           }
           CCTK_REAL invw = 1./sqrt( 1. + psi4*SQR(uIso) );

           rho[idx] = rhoL;
           eps[idx] = eos_K * pow( rhoL, eos_gamma-1. ) / ( eos_gamma - 1.);
           velx[idx] = -uIso * (xL/rIso) * invw;
           vely[idx] = -uIso * (yL/rIso) * invw;
           velz[idx] = -uIso * (zL/rIso) * invw;

        }
     }
  }

}

static CCTK_INT calc_BondiSolution( CCTK_REAL r, CCTK_REAL * rho, CCTK_REAL * u, 
   const CCTK_REAL rSonic, const CCTK_REAL rhoSonic, 
   const CCTK_REAL M, const CCTK_REAL Mdot, const CCTK_REAL eos_K, 
   const CCTK_REAL Qdot ) {

   DECLARE_CCTK_PARAMETERS;

   /* Store initial guesses ... */
   CCTK_REAL rho_0, u_0;
   rho_0 = *rho;
   u_0 = *u;

   /* Newton-Raphson Parameters */
   const CCTK_REAL max_its = 50;
   const CCTK_REAL newt_tol = 1.e-16;
   const CCTK_REAL sol_tol =  5.e-16; //0.7*newt_tol;
   const CCTK_INT  extra_iterations = 2;
   const CCTK_INT  num_rounds_max = 1;
   const CCTK_REAL reset_rho = 1e-5*rho_min;

   /* Case where we had a bad initial rho. */
   if ( rho_0 < 0. ) {
      if ( r < 1.1*rSonic && r > 0.9*rSonic ) {
         r = rSonic;
      } else {
         CCTK_REAL u_r;
         if ( r<rSonic ) {
            u_r = 1./sqrt(r);
         } else {
            u_r = 0.5*pow(r, -1.5);
         }
         rho_0 = Mdot/(4.*M_PI*SQR(r)*u_r);
      }
   }

   /* ... using Newton-Raphson, ... */
   /* Note: there are 2 branches of solutions we might end up on. */
   CCTK_REAL rho_old;
   CCTK_REAL rho_tmp = rho_0;
   CCTK_REAL u_tmp = u_0;
   CCTK_INT  itn_count = 0;
   CCTK_INT  num_rounds = 0;
   CCTK_INT  itn_count_extra = 0;
   CCTK_INT  satisfied = 0;
   CCTK_INT  keep_going = 1;
   CCTK_INT  doing_extra = 0;
   CCTK_REAL err = 1e3;
   CCTK_REAL QdotL, errQdot, MdotL, errMdot;


   while ( keep_going==1 ) {

      /* Travel along Mdot=const */
      CCTK_REAL h, dhdrho, dudrho, uterm, resid;
      dhdrho = eos_K * eos_gamma * pow( rho_tmp, eos_gamma-2. );
      //h = 1. + eos_K * eos_gamma * pow( rho_tmp, eos_gamma-1.)/(eos_gamma - 1.);
      h = 1. + dhdrho*rho_tmp/(eos_gamma - 1.);
      u_tmp = Mdot / ( 4. * M_PI * SQR(r) * rho_tmp );
      dudrho = -u_tmp/rho_tmp;
      uterm = 1. - 2.*M/r + SQR(u_tmp);

      CCTK_REAL f, df, drho;
      f = SQR(h) * ( 1. - 2.*M/r + SQR(u_tmp) ) - Qdot;
      df = 2.*h*( dhdrho*uterm + h*u_tmp*dudrho );
      drho = -f/df; 


      /* Save & Update rho */
      rho_old = rho_tmp;
      rho_tmp += drho;

      /* Error: Absolute error in Qdot(rho,u) */
      QdotL = SQR(h)*(1.-2.*M/r + SQR(u_tmp)); 
      errQdot = (QdotL - Qdot)/Qdot;
      MdotL = 4.*M_PI*SQR(r)*rho_tmp*u_tmp;
      errMdot = (MdotL-Mdot)/Mdot;

      /* Error and convergence: relative change in rho over iteration */
      if ( rho_tmp <= 0 ) {
         err = fabs(drho);
         rho_tmp = reset_rho;
         //if ( debug ) CCTK_VInfo(CCTK_THORNSTRING, "   Reset rho at iteration %d-%d.",num_rounds, itn_count);
      } else {
         err = fabs(drho/rho_tmp);
      }

      /* Relative error satisfied. Extra iterations & extra rounds. */
      if ( ( fabs(err) < newt_tol || fabs(f) < sol_tol ) && doing_extra==0 && extra_iterations>0  ) {
         satisfied = 1;
         doing_extra = 1;
      }
      if ( doing_extra == 1 ) itn_count_extra++;

      /* Reset/check keep_going if we're to head into (are in) extra iterations */
      if ( (satisfied==1 && doing_extra==0) || ( itn_count_extra > extra_iterations) ) { 
         keep_going=0;
      }

      /* Start a new round */
      if ( itn_count+1>=max_its && fabs(err) > newt_tol && num_rounds<num_rounds_max ) {
         num_rounds++;
         itn_count=-1;
         rho_tmp=1.05*rho_tmp; /* tweak rho to hopefully get better convergence */
      }

      /* Generic fall out if itn & group maxes are met */
      if ( (fabs(err) > newt_tol || fabs(f) > sol_tol) && num_rounds+1>=num_rounds_max && itn_count+1>=max_its ) {
         keep_going=0;
      }

      itn_count++;

      /* Very verbose */
      /*if ( debug ) {
         CCTK_VInfo(CCTK_THORNSTRING, "     Itn (%d-%d): (rho,u)=(%g,%g)->(%g,%g), drho=%g, (f,df)=(%g,%g), err=%g; delta(Mdot,Qdot)=(%g,%g)",
             num_rounds, itn_count, rho_old, u_tmp, rho_tmp, Mdot/(4.*M_PI*SQR(r)*rho_tmp), drho, f, df, err, errMdot, errQdot);
      }*/
   }

   if ( satisfied == 0 ) {

      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
          "Failed to find Bondi solution for r=%g, EoS=(%g,%g), Mdot=%g, Qdot=%g, rSonic=%g, (rho,u)=(%g,%g)->(%g,%g) with err=%g; delta(Mdot,Qdot)=(%g,%g).",
          r, eos_K, eos_gamma, Mdot, Qdot, rSonic, rho_0, u_0, rho_tmp, u_tmp, err, errMdot, errQdot); 

   }
     /* Replace initial guess with new solution */
     //*rho = rho_tmp;
     //*u   = Mdot / ( 4. * M_PI * SQR(r) * rho_tmp );
   //}

   *rho = rho_tmp;
   *u   = Mdot / ( 4. * M_PI * SQR(r) * rho_tmp );
   return satisfied;

}

void SphericalInfall_IsoBondi_WriteSolns(CCTK_ARGUMENTS)
{

   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   char out_format[5];
   sprintf(out_format, ".19g");

   for ( int nsol=0; nsol<bondi_num_solutions; nsol++ ) {

       size_t length_written; 
       FILE *file;
       char *filename;

       assert (r_bondi[bondi_maxnpts*nsol]);

       Util_asprintf ( &filename, "%s/Bondi_Solution_%d.asc", out_dir, nsol );
       assert(filename);

       file = fopen (filename,"w");
       if (!file) {
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
              "Could not open bondi solution output file '%s'", filename); 
          free(filename);
          continue;
       }
       
       /* Header information */
       fprintf(file, "# Bondi Solution %d\n", nsol);
       fprintf(file, "#   M = %g at (%g,%g,%g)\n", 
               bondi_central_mass[nsol], bondi_bh_x0[nsol], bondi_bh_y0[nsol],
               bondi_bh_z0[nsol]); 
       fprintf(file, "#   Mdot = %g, r_sonic = %g\n", bondi_mdot[nsol], bondi_sonicpt_radius[nsol] ); 
       fprintf(file, "#   %d points, r = [ %g : %g ]\n", bondi_npts[nsol], bondi_rmin[nsol], bondi_rmax[nsol]); 
       fprintf(file, "# Data headers:\n");
       fprintf(file, "# 1:r(Iso) 2:r(Sch) 3:rho 4:u_Sch 5:v_Sch 6:v_Iso 7:delta(Mdot) 8:delta(Qdot)\n");
      

       /* Print data */
       char format_str[2048];
       size_t len_written = snprintf (format_str, sizeof(format_str)/sizeof(format_str[0]),
                               "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", out_format, out_format, out_format,
                               out_format, out_format, out_format, out_format, out_format);
       assert( len_written < sizeof(format_str)/sizeof(format_str[0])); 
       CCTK_REAL M = bondi_central_mass[nsol];
       for ( int i=0; i<bondi_npts[nsol]; i++ ) {

           /* Solution */
           CCTK_REAL rSch = r_bondi[i];
           CCTK_REAL uSch = u_bondi[i];

           /* Derive other quantities for output */
           CCTK_REAL rIso, psi, vSch, vIso, uIso, invwIso;
           if ( rSch >=2.*M ) {
              rIso = 0.5*( sqrt( SQR(rSch) - 2.*M*rSch) + rSch - M );
              psi = 1. + M/(2.*rIso);
           } else {
              /* FIX */
              rIso = 0;
              psi = 1.;
           }
           vSch = -uSch / sqrt( 1. - 2.*M/rSch + SQR(uSch) );
           uIso = uSch / ( 1. - SQR( M/(2.*rIso) ) );
           invwIso = 1./sqrt( 1. + QUAD(psi)*SQR(uSch) );
           vIso = uIso * invwIso; 

           /* Fractional error, Mdot */
           CCTK_REAL errMdot = (4.*M_PI*SQR(rSch)*rho_bondi[i]*uSch - bondi_mdot[nsol])/bondi_mdot[nsol];

           /* Fractional error, Qdot */
           CCTK_REAL errQdot, Qdot, QdotL, eos_K, aSonic2, hSonic, rhoSonic;
           aSonic2 = bondi_central_mass[nsol]/ ( 2.*bondi_sonicpt_radius[nsol] - 3.*bondi_central_mass[nsol] );
           hSonic = 1./( 1. - aSonic2/(eos_gamma-1.));
           rhoSonic = bondi_mdot[nsol] / ( 4. * M_PI * sqrt(bondi_central_mass[nsol]/(2.*bondi_sonicpt_radius[nsol])) * SQR(bondi_sonicpt_radius[nsol]) );
           Qdot = (1. - 1.5*bondi_central_mass[nsol]/bondi_sonicpt_radius[nsol]) * SQR(hSonic);
           eos_K = hSonic * aSonic2 * pow( rhoSonic, 1.-eos_gamma ) / eos_gamma;
           QdotL = SQR( 1. + eos_gamma*eos_K*pow(rho_bondi[i],eos_gamma-1.)/(eos_gamma-1.))*( 1. - 2.*bondi_central_mass[nsol]/rSch + SQR(uSch) );
           errQdot = ( QdotL - Qdot ) / Qdot;
 

           /* Final output */ 
           fprintf(file, format_str, rIso, rSch, rho_bondi[i], uSch, vSch, vIso, errMdot, errQdot ); 

       }
       fprintf(file,"\n");

       fclose(file);
       free(filename);

   }

   /* DEBUG! */
   //CCTK_WARN(0,"Check Bondi Solution output before we proceed further.");

}
