/* Handwritten solve of con 2 prim 
 * http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "EOS_Base.h"
#include <stdio.h>
#include <math.h>

#define JUNK -1e3

/***************************/
/* Primary solver routines */
/***************************/

int EOS_solveHDsys (int eos_handle,
    bool is_polytype, double min_rho, double pmin,
    double udens, double utau, double usx, double usy, double usz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *rho, double *eps, double *ux, double *uy, double *uz, double *soln_quality)
{

   DECLARE_CCTK_PARAMETERS;
   int status;

   if ( verbose > 6 ) {
	CCTK_VInfo(CCTK_THORNSTRING,"Entered %sEOS_solvesys with the following quantities ...", is_polytype ? "polytrope" : "general");
	CCTK_VInfo(CCTK_THORNSTRING,"    g_ij = (%g,%g,%g,%g,%g,%g)",gxx,gxy,gxz,gyy,gyz,gzz);
	CCTK_VInfo(CCTK_THORNSTRING,"    g^ij = (%g,%g,%g,%g,%g,%g)",guxx,guxy,guxz,guyy,guyz,guzz);
   }

   // CCTK_WARN(0,"HD via MHD with HD recovery is not functional at the moment. Please debug.");

   CCTK_REAL rel_err = 999;
   CCTK_REAL abs_err = 999;
   CCTK_REAL atmo_reset = soln_quality[2];

   /* Constants in our iterations */
   CCTK_REAL uS2, eps_guess;
   uS2 = SQR(usx)*guxx + SQR(usy)*guyy + SQR(usz)*guzz
	+ 2.*usx*usy*guxy + 2.*usx*usz*guxz + 2.*usy*usz*guyz;

   if ( ! is_polytype ) {

      /* Big loop -- solve for f=P-P */
      eps_guess = *eps;
      CCTK_REAL p_old = fmax(( 1.+press_tolerance)*pmin, EOS_Pressure(eos_handle,*rho,eps_guess));
      CCTK_REAL p_low = fmax(pmin, sqrt(uS2) - utau - udens );

      if ( SQR(utau + udens + p_old) - uS2 <= 0 ) { // Check if initial guess is physical
         if (reset_pressure != 0) {
            p_old = pmin + ( sqrt(uS2) - utau - udens );
         } else {
            *eps = -1.; // Poison
            return 1; 
         }
      }

      /* Find rho, W, epsilon given conservatives and pressure guess */
      CCTK_REAL it_rho, w_lorentz, epsilon;
      it_rho = udens * sqrt( SQR(utau+udens+p_old) - uS2)/(utau+p_old+udens);
      w_lorentz = 1./ sqrt( 1. - uS2/SQR(utau+p_old+udens));
      epsilon = ( sqrt( SQR(utau+udens+p_old) - uS2) - p_old*w_lorentz - udens)/udens;

      /* Handwritten iterative solver */
      CCTK_REAL dp_drho, dp_deps, tmp_tpds, drho_dp, deps_dp, df;
      CCTK_REAL f = p_old - EOS_Pressure(eos_handle, it_rho, epsilon);
      int count=0;
      CCTK_REAL p_new = p_old;

      if( udens < min_rho*(1.+atmo_tolerance) ) { // Skip if we're going to be resetting to atmo anyways 
        // Actually solve
        while ( ( ( rel_err > solver_epsrel ) && ( abs_err > solver_epsabs ) ) || count < 1 )
        {
            count = count + 1.;
            if ( count > countmax ) {
              status=1;
              *eps = -1;
              break;
            }

            dp_drho = EOS_DPressByDRho(eos_handle, it_rho, epsilon);
            dp_deps = EOS_DPressByDEps(eos_handle, it_rho, epsilon);
            tmp_tpds = SQR(utau+udens+p_new)-uS2;
            drho_dp = udens*uS2/( sqrt(tmp_tpds)*SQR(utau+udens+p_new) );
            deps_dp = p_new * uS2 / ( it_rho * (udens+utau+p_new)*tmp_tpds );
            
            df = 1. - dp_drho*drho_dp - dp_deps*deps_dp;

            p_old = p_new;
            p_new = p_old - f/df; /* Newton-Raphson standard */

            /* If p would go negative, use fake bisection. */
            tmp_tpds = SQR( utau+udens+p_new) - uS2;
            if ( p_new <= p_low || tmp_tpds <= 0 ) { 
               p_new = 0.5*(p_old+p_low); 
               tmp_tpds = SQR( utau+udens+p_new) - uS2;
               if ( p_new <= p_low || tmp_tpds < 0 ) { // Fake bisection failed.
                  if ( verbose > 3 && p_new < p_low )
                      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                                 "Con2Prim: Pressure invalid even after pressure reset (%g<%g) tmp=%g, f=%g,df=%g",
                                 p_new,p_low,tmp_tpds,f,df);
                  p_new = p_old; // This breaks loop.
                  tmp_tpds = SQR( utau+udens+p_new) - uS2;
                  if ( SQR(utau+udens+p_new)-uS2 <= 0 )
                     CCTK_WARN(0,"Con2Prim failed completely.");
                  // Currently omit extra verbosity 
               }
            }

            /* New p safe, we recalculate f and company */
            it_rho = udens*sqrt(tmp_tpds)/(utau+udens+p_new);
            w_lorentz = 1./sqrt( 1. - uS2/SQR(utau+udens+p_new));
            epsilon = ( sqrt(tmp_tpds) - p_new*w_lorentz - udens)/udens;

            f = p_new - EOS_Pressure(eos_handle, it_rho, epsilon);

            abs_err = fabs(p_new-p_old);
            rel_err = abs_err/fabs(p_new);

        } // end loop
        abs_err = fabs(p_new-p_old);
        rel_err = abs_err/fabs(p_new);

        // We skip root polishing since we don't tend to use that option.
        status = ( rel_err > solver_epsrel && abs_err > solver_epsabs ) ? 1 : 0; // success=0

      } else {
        // Atmo is atmo
        status = 0;
        abs_err = 0;
        rel_err = 0;
      } // End non-atmo specific 

      soln_quality[0] = rel_err;
      soln_quality[1] = abs_err;
      

      if( status==0 ) {
        if( ( udens < min_rho*(1.+atmo_tolerance) || it_rho < min_rho*(1.+atmo_tolerance) ) ) {
           *rho = min_rho;
           *ux = 0.;
           *uy = 0.;
           *uz = 0.;
           *eps = EOS_SpecificIntEnergy( eos_handle, min_rho, p_new );
           atmo_reset += 1.;
           /* Reset soln_quality */
        } else {
          // We have rho, w, epsilon.
          *rho = it_rho;
          CCTK_REAL rho_h_W = ( it_rho*(1. + epsilon) + p_new )*w_lorentz;
          /* u_i = S_i / (rho*h*W) */
          *ux= usx/rho_h_W;
          *uy= usy/rho_h_W;
          *uz= usz/rho_h_W;
          *eps = epsilon;
        }
      }
      soln_quality[2] += atmo_reset;
      return status;

 
   } else { // if polytype, solve for f=D-rho W instead //

      CCTK_REAL epsilon = JUNK;
      CCTK_REAL rho_old, press, enthalpy, w_lorentz, f, df;

      /* Initial state from rho, which should be last primitive */
      rho_old = fmax( min_rho, *rho );
      press = EOS_Pressure(eos_handle, rho_old, JUNK);
      epsilon = EOS_SpecificIntEnergy(eos_handle, rho_old, press);
      enthalpy = 1. + EOS_SpecificIntEnergy(eos_handle, rho_old, press) + press/rho_old; 
      w_lorentz = sqrt( 1. + uS2/SQR(udens*enthalpy) );

      f = rho_old * w_lorentz - udens;

      int count=0, status=0;
      CCTK_REAL rho_new = rho_old; 
      while ( ((rel_err > solver_epsrel) && 
                (abs_err > solver_epsabs)) || count < 1 ) {

          count = count + 1;
          if ( count > countmax ) {
             status=0;
             break;
          }

          CCTK_REAL denthalpy = EOS_DPressByDRho( eos_handle, rho_new, epsilon ) / rho_new;
          df = w_lorentz - rho_new*uS2*denthalpy / ( w_lorentz*SQR(udens)*SQR(enthalpy)*enthalpy );

          rho_old = rho_new;
          rho_new = rho_old - f/df;

          /* Impose positivity, approach low-density more slowly. */
          if ( rho_new <= 0 )
             rho_new = rho_old / 2.;

          /* Recalculate values */
          epsilon = EOS_SpecificIntEnergy( eos_handle, rho_new, press );
          press = EOS_Pressure(eos_handle, rho_new, epsilon );
          enthalpy = 1. + epsilon + EOS_Pressure(eos_handle, rho_new, epsilon )/rho_new;
          w_lorentz = sqrt( 1. + uS2/SQR(udens*enthalpy) );

          f = rho_new*w_lorentz - udens;
          abs_err = fabs(rho_new-rho_old);
          rel_err = abs_err/fabs(rho_new);

      } // End loop
      abs_err = fabs(rho_new-rho_old);
      rel_err = abs_err/fabs(rho_new);

      status = ( rel_err > solver_epsrel && abs_err > solver_epsabs ) ? 1 : 0; // success=0
      soln_quality[0] = rel_err;
      soln_quality[1] = abs_err;

      if ( status == 0 ) {
         if( udens < min_rho*(1.+atmo_tolerance) || rho_new < min_rho*(1.+atmo_tolerance) ) {
            *rho = min_rho;
            *ux = 0.;
            *uy = 0.;
            *uz = 0.;
            atmo_reset = atmo_reset + 1.;
         } else {
            CCTK_REAL rho_h_W = w_lorentz*( rho_new*(1. + epsilon) + press );
            *rho = rho_new;
            *ux= usx/rho_h_W;
            *uy= usy/rho_h_W;
            *uz= usz/rho_h_W;
        }
        soln_quality[2] = atmo_reset;
        return status;
      }

   }

}
