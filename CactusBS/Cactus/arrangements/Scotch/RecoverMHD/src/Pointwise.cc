/* Functions for converting conservatives to primitives
   on a pointwise basis */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "EOS_Base.h"
#include "RecoverMHD.h"
#include "cctk_WarnLevel.h"
#include <assert.h>
#define MINB2 -1 //1e-30
#define FUDGEEPS 1e-10
#define JUNK -1e3

// Include non-gsl solver code for speed 
#include "Solve_HDSys.cc"
#include "Solve_v2Z.cc"
#include "Solve_W.cc"

extern "C" 
CCTK_INT RecoverMHD_Con2Prim_pointwise( CCTK_INT handle, CCTK_INT phandle, CCTK_INT polytype, CCTK_INT iteration, CCTK_INT stype,
	CCTK_INT sys_eqns, CCTK_REAL *dens, CCTK_REAL *tau, CCTK_REAL *sx, CCTK_REAL *sy, CCTK_REAL *sz, CCTK_REAL *Bnx, CCTK_REAL *Bny, CCTK_REAL *Bnz,
	CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz, CCTK_REAL *Bvecx, CCTK_REAL *Bvecy, CCTK_REAL *Bvecz,
	CCTK_REAL *w_lorentz, CCTK_REAL *press, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
	CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
	CCTK_INT whisky_reflevel, CCTK_REAL carpetweight, 
	CCTK_INT do_impose_s2_limit, CCTK_REAL rho_min, CCTK_REAL tau_min )
{
   DECLARE_CCTK_PARAMETERS;

   /* GSL Solver type setup */
   typedef enum GSL_SOLVER_TYPE gsolve_type;
   gsolve_type gstype;
   gstype = GST_Newton; // HARD-CODE Newton-Raphson
   /*switch ( stype ) {
	case 0:
		gstype = GST_Newton;
		break;
	case 1:
		gstype = GST_gNewton;
		break;
	case 2:
		gstype = GST_HybridSJ;
		break;
	case 3:
		gstype = GST_HybridJ;
		break;
	case 4:
		gstype = GST_HybridS;
		break;
	case 5:
		gstype = GST_Hybrid;
		break;
	case 6:
		gstype = GST_dNewton;
		break;
	case 7:
		gstype = GST_Broyden;
		break;
        default:
                CCTK_WARN(0,"Invalid GSL Multiroot solver handle in Con2PrimMHD call.");
                gstype = (gsolve_type)0; // make gcc happy //
                break;
   } */

   /* Save initial values for verbosity */
   CCTK_REAL dens_in, tau_in, sx_in, sy_in, sz_in, Bnx_in, Bny_in, Bnz_in;
   CCTK_REAL rho_in, eps_in, velx_in, vely_in, velz_in, Bvx_in, Bvy_in, Bvz_in;
   dens_in = *dens;
   sx_in = *sx;
   sy_in = *sy;
   sz_in = *sz;
   Bnx_in = *Bnx;
   Bny_in = *Bny;
   Bnz_in = *Bnz;
   rho_in = *rho;
   eps_in = *eps;
   velx_in = *velx;
   vely_in = *vely;
   velz_in = *velz;
   Bvx_in = *Bvecx;
   Bvy_in = *Bvecy;
   Bvz_in = *Bvecz;
   tau_in = *tau;

   /* Toggle whether the conservatives have changed at all */ 
   bool cons_have_changed=0; 

   /* Determinant and Inverse metric */
   CCTK_REAL ugxx, ugxy, ugxz, ugyy, ugyz, ugzz, det, rootdet, inv_sqrt_det, psiL;
   det = spatialdeterminant( gxx, gxy, gxz, gyy, gyz, gzz );
   invertspatialdet( gxx, gxy, gxz, gyy, gyz, gzz, det,
		ugxx, ugxy, ugxz, ugyy, ugyz, ugzz );
   if ( carpetweight == 1 && (std::isnan(det) || det < 0) ) {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Con2PrimMHD: Metric determinant in iteration %d "
		"at (%g,%g,%g) on level %d is: %g from data g=(%g,%g,%g,%g,%g,%g)",
		iteration, x,y,z,whisky_reflevel,det,gxx,gxy,gxz,gyy,gyz,gzz);
    }
    rootdet = sqrt(det);
    inv_sqrt_det = 1/rootdet;
    psiL=pow( rootdet, 1./6. );


    //if ( verbose > 8 ) {
    if ( verbose > 8 || fabs(*Bnx) > 1e20 ) {
	CCTK_INFO("Con2PrimMHD has been called with the following parameters:");
	CCTK_VInfo(CCTK_THORNSTRING,"   Metric: g_ij=(%g,%g,%g,%g,%g,%g) with psiL=%g)",gxx,gxy,gxz,gyy,gyz,gzz,psiL);
	CCTK_VInfo(CCTK_THORNSTRING,"           ( detg=%g, g^ij=(%g,%g,%g,%g,%g,%g)) as just calculated.",det,ugxx,ugxy,ugxz,ugyy,ugyz,ugzz);
	CCTK_VInfo(CCTK_THORNSTRING,"   Hydro : dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g)",*dens,*sx,*sy,*sz,*Bnx,*Bny,*Bnz);
	CCTK_VInfo(CCTK_THORNSTRING,"           tau=%g",*tau);
	
    }

    /* Extra limits */    
    CCTK_REAL pmin = EOS_Pressure( phandle, rho_min, JUNK );
    CCTK_REAL epsmin = EOS_SpecificIntEnergy( phandle, rho_min, pmin );
    CCTK_REAL atmo_threshold = rho_min*(1.+atmo_tolerance);

    /* Initialize quality controls */
    bool generalEOS_unsat=0;
    int fullstatus=0;
    int ierr=0;
    CCTK_REAL soln_quality[3];
        soln_quality[0] = 1000;
        soln_quality[1] = 1000;
        soln_quality[2] = 0.;

    /***** Recover B-field *****/
    /********* Trivial *********/
    *Bvecx = *Bnx * inv_sqrt_det;
    *Bvecy = *Bny * inv_sqrt_det;
    *Bvecz = *Bnz * inv_sqrt_det;

    /* Useful magnitudes */
    CCTK_REAL S2 = SQR(*sx)*ugxx + SQR(*sy)*ugyy + SQR(*sz)*ugzz
	+ 2.*(*sx)*(*sy)*ugxy + 2.*(*sx)*(*sz)*ugxz + 2.*(*sy)*(*sz)*ugyz;
    CCTK_REAL B2 = SQR(*Bvecx)*gxx + SQR(*Bvecy)*gyy + SQR(*Bvecz)*gzz 
        + 2*( gxy*(*Bvecx)*(*Bvecy) + gxz*(*Bvecx)*(*Bvecz) + gyz*(*Bvecy)*(*Bvecz) );
    CCTK_REAL Bn2 = SQR(*Bnx)*gxx + SQR(*Bny)*gyy + SQR(*Bnz)*gzz
        + 2.*gxy*(*Bnx)*(*Bny) + 2.*gxz*(*Bnx)*(*Bnz) + 2.*gyz*(*Bny)*(*Bnz);
    if ( B2 < 0 ) 
       CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
       "B^2 < 0 in recovery at (%g,%g,%g) with detg=%g",
       x,y,z,det);

    /*********************************************/
    /*** Conservative Checks                   ***/
    /*** Modified from Etienne arXiv:1112.0568 ***/
    /*** Changes are considered when           ***/
    /***   dens <= 0                           ***/
    /***   tau < (1/2)B^2                      ***/
    /***   other checks w/i large psi regions  ***/
    /*********************************************/

    if ( *dens <= 0 ) {

       /* Already set to atmosphere here! */
       *rho = rho_min;
       *press = pmin;
       *eps = epsmin; 
       *velx  = 0.; 
       *vely  = 0.; 
       *velz  = 0.;

       Prim2ConMHDPoly( phandle, gxx, gxy, gxz, gyy, gyz, gzz, det,
                        dens, sx, sy, sz, tau, *rho, *velx, *vely, *velz,
                        eps, press, w_lorentz, Bnx, Bny, Bnz, 
                        *Bvecx, *Bvecy, *Bvecz );
       /* No use going on! */
       return fullstatus;
 
       cons_have_changed = 1;
       soln_quality[2] = soln_quality[2] + 1.;

    }

    /* Check tau >= (1/2)*Bn2/sqrt(det) (Remember, we have engulfed the 1/4pi in the def of B) */
    /* As in, tau = tau_g + tau_mag. */
    if ( tau != NULL ) {
       if ( *tau < 0.5*inv_sqrt_det*Bn2 ) {
          CCTK_REAL old_tau = *tau;
          *tau = 0.5*rootdet*B2 + tau_min;
          if ( verbose > 2 ) {
             CCTK_REAL v2_in = gxx*SQR(velx_in)+gyy*SQR(vely_in)+gzz*SQR(velz_in)
                        + 2*gxy*velx_in*vely_in + 2*gxz*velx_in*velz_in + 2*gyz*vely_in*velz_in;
             CCTK_REAL w_in = 1./sqrt( 1. - v2_in );
             if ( *dens > rootdet*rho_min*w_in ) {
                CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, 
                        "Reset tau from %g to %g=%15.7g*old_tau at r=(%g,%g,%g)=%g!",
                        old_tau, *tau, *tau/old_tau, x,y,z, sqrt(SQR(x)+SQR(y)+SQR(z)));
             } 
          } 
          cons_have_changed = 1;
       } 
    }

    /* Apply cap to magnetic field first, since this also changes the momentum considerably. */
    if ( do_impose_s2_limit ) { /* This toggle should only exist within a BH's horizon */

       /* Apply a dust limit if necessary to cap unphysically large momenta */
       if (!polytype) {
          CCTK_REAL maxS2 = recovery_flimit * (*tau) * (*tau + 2.*(*dens));
          if ( S2 > maxS2 ) {
             // limit appears twice
             CCTK_REAL flimit = sqrt(recovery_flimit * maxS2/S2);
             if ( std::isnan(flimit) && carpetweight == 1 ) {
                CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Scale factor is NaN! (tau, dens, S2, maxS2)=(%g,%g,%g,%g) at (%g,%g,%g)",
                           *tau, *dens, S2, maxS2, x, y, z);
             }
             //CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Trying to implement dust limit in MHD at (%g,%g,%g). Please decide on a correction!",x,y,z);
             *sx = flimit * *sx;
             *sy = flimit * *sy;
             *sz = flimit * *sz;
             S2 = SQR(*sx)*ugxx + SQR(*sy)*ugyy + SQR(*sz)*ugzz
                  + 2.*(*sx)*(*sy)*ugxy + 2.*(*sx)*(*sz)*ugxz + 2.*(*sy)*(*sz)*ugyz;
             cons_have_changed = 1;
          }
       }

       /* Impose an extra |B| cap */ 
       if ( impose_b2_limit ) {
          if ( Bn2 > SQR(maxB) ) {
             // Apply limit to both primitive and conservative, and 
             CCTK_REAL flimit = maxB/sqrt(Bn2);

             if (verbose>1) {
                CCTK_VInfo(CCTK_THORNSTRING,"Rescaled Bn from (%g,%g,%g) to (%g,%g,%g)",
                           *Bnx, *Bny, *Bnz, flimit * *Bnx, flimit * *Bny, flimit * *Bnz );
             }
             *Bnx = flimit * *Bnx;
             *Bny = flimit * *Bny;
             *Bnz = flimit * *Bnz;

             cons_have_changed = 1;

             /* Update Bvec, Bn2 and B2 */
             *Bvecx = *Bnx * inv_sqrt_det;
             *Bvecy = *Bny * inv_sqrt_det;
             *Bvecz = *Bnz * inv_sqrt_det;
             B2 = SQR(*Bvecx)*gxx + SQR(*Bvecy)*gyy + SQR(*Bvecz)*gzz 
                 + 2*( gxy*(*Bvecx)*(*Bvecy) + gxz*(*Bvecx)*(*Bvecz) 
                 + gyz*(*Bvecy)*(*Bvecz) );
             Bn2 = SQR(*Bnx)*gxx + SQR(*Bny)*gyy + SQR(*Bnz)*gzz
                 + 2.*gxy*(*Bnx)*(*Bny) + 2.*gxz*(*Bnx)*(*Bnz) 
                 + 2.*gyz*(*Bny)*(*Bnz);

             // Need to keep s and tau consistent. Use old w, this will underestimate the contribution of course.
             //CCTK_REAL inv_w_in = 1. - ( SQR(velx_in) + SQR(vely_in) + SQR(velz_in) );
             //*sx = *sx - (1.-flimit) * ( ) * inv_w_in

          } 
       }
    }

    /* Are we deep within a BH? */
    if ( (rootdet >= pow(deep_BH_psi,6)) && B2>0 ) {

       CCTK_REAL SdotBhat = ( (*sx)*(*Bvecx) + (*sy)*(*Bvecy) + (*sz)*(*Bvecz) )/B2;
       CCTK_REAL Wm0 = inv_sqrt_det * sqrt( SQR(SdotBhat) + SQR(*dens) );
       CCTK_REAL Sm02 = ( SQR(Wm0)*S2 + SQR( SdotBhat )*( B2 + 2*Wm0 ) )/SQR( Wm0 + B2 );
       CCTK_REAL Wm = inv_sqrt_det*sqrt( Sm02 + SQR(*dens) );

       CCTK_REAL Sm = sqrt( (SQR(Wm)*S2 + SQR(SdotBhat)*(B2+2.*Wm))/SQR(Wm+B2) );
       CCTK_REAL taum = (*tau) - 0.5*rootdet*B2 - 0.5*inv_sqrt_det*(B2*S2 - SQR(SdotBhat))/SQR(Wm+B2);

       if ( taum < tau_min ) {
          CCTK_REAL old_tau = *tau;
          *tau = tau_min + 0.5*rootdet*B2 + 0.5*inv_sqrt_det*( B2*S2 - SQR(SdotBhat) )/SQR(Wm+B2);
          //CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, 
          //           "Modified tau from %g to %g at r=(%g,%g,%g) due to tau_m calculation!",
          //           old_tau, *tau, x,y,z);
          cons_have_changed = 1;
       }

       CCTK_REAL maxS2 = recovery_flimit * taum*(taum + 2.*(*dens));
       if ( S2 > maxS2 ) {
          //CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, 
          //           "Modifying S2 from %g to %g at r=(%g,%g,%g) due to tau_m calculation!",
          //           S2, maxS2, x,y,z);
          CCTK_REAL limitfac = sqrt( maxS2/S2 );
          *sx = limitfac*(*sx);
          *sy = limitfac*(*sy);
          *sz = limitfac*(*sz);
          S2 = SQR(*sx)*ugxx + SQR(*sy)*ugyy + SQR(*sz)*ugzz
               + 2.*(*sx)*(*sy)*ugxy + 2.*(*sx)*(*sz)*ugxz + 2.*(*sy)*(*sz)*ugyz;
          cons_have_changed = 1;
       }

    }

    /* DONE possibly manipulating incoming conservatives */ 
    if ( cons_have_changed && verbose > 2 && carpetweight > 0.5 ) {
       CCTK_VInfo(CCTK_THORNSTRING,"Conservatives Modification (polytype=%d) at x=(%g,%g,%g), psi=%g, rl=%d, weight=%g:",
                  polytype,x,y,z,pow(psiL,0.25),whisky_reflevel,carpetweight);
       CCTK_VInfo(CCTK_THORNSTRING,"   dens  = %g --> %g", dens_in, *dens);
       CCTK_VInfo(CCTK_THORNSTRING,"   tau   = %g --> %g", tau_in, *tau);
       CCTK_VInfo(CCTK_THORNSTRING,"   sx    = %g --> %g", sx_in, *sx);
       CCTK_VInfo(CCTK_THORNSTRING,"   sy    = %g --> %g", sy_in, *sy);
       CCTK_VInfo(CCTK_THORNSTRING,"   sz    = %g --> %g", sz_in, *sz);
       CCTK_VInfo(CCTK_THORNSTRING,"   Bnx   = %g --> %g", Bnx_in, *Bnx);
       CCTK_VInfo(CCTK_THORNSTRING,"   Bny   = %g --> %g", Bny_in, *Bny);
       CCTK_VInfo(CCTK_THORNSTRING,"   Bnz   = %g --> %g", Bnz_in, *Bnz);
    }

    /*************************************
     *** Finish set-up for prim solver ***
     *************************************/

    /* Undensitize conservative variables */
    CCTK_REAL udens, usx, usy, usz, s2, utau;
    udens = *dens *inv_sqrt_det;
    usx = *sx * inv_sqrt_det;
    usy = *sy * inv_sqrt_det;
    usz = *sz * inv_sqrt_det;
    s2 = SQR(usx)*ugxx + SQR(usy)*ugyy + SQR(usz)*ugzz
	+ 2.*usx*usy*ugxy + 2.*usx*usz*ugxz + 2.*usy*usz*ugyz;
    utau = *tau * inv_sqrt_det;

    /* Initial Plasma beta parameter */
    CCTK_REAL vlowx, vlowy, vlowz;
    vlowx = gxx*(*velx) + gxy*(*vely) + gxz*(*velz);
    vlowy = gxy*(*velx) + gyy*(*vely) + gyz*(*velz);
    vlowz = gxz*(*velx) + gyz*(*vely) + gzz*(*velz);
    CCTK_REAL b2_comoving_in = B2/SQR(*w_lorentz) + SQR(vlowx*(*Bvecx) + vlowy*(*Bvecy) + vlowz*(*Bvecz));
    CCTK_REAL betaPlasma_in = 0.5*b2_comoving_in/(*press);

    /* Go through the solved system type only once per point (in case of fallback) */
    /* Speed things up and move this out of this routine in the future? */
    enum SYSTEM_EQN_TYPE{ HD_ONLY=0, MHD_5D_UiEps=1, MHD_2D_WP=2, MHD_1D_hW2=3, MHD_2D_v2Z=4 } solve_system_method;
    if ( B2 < MINB2 ) {
       solve_system_method = HD_ONLY;
    } else { 
       solve_system_method = (SYSTEM_EQN_TYPE) sys_eqns;
    }

    /* Guess conversion by old primitives for initial solver state */
    /* Recover hydro variables */
    CCTK_REAL eps_guess, ux_guess, uy_guess, uz_guess; /* u[xyz] = u_i */
    CCTK_REAL rho_new, eps_new, ux_new, uy_new, uz_new;
    if ( polytype ) {
	eps_guess = EOS_SpecificIntEnergy(handle,*rho,*press);
    } else {
	eps_guess = *eps;
    }
    ux_guess = vlowx * (*w_lorentz);
    uy_guess = vlowy * (*w_lorentz);
    uz_guess = vlowz * (*w_lorentz);

    ux_new = ux_guess;
    uy_new = uy_guess;
    uz_new = uz_guess;

    rho_new = rho_in; 
    eps_new = eps_guess;

    /* Don't even bother with the iterative method if we're going to apply atmo anyways */
    if ( udens < atmo_threshold ) {

       *rho = rho_min;
       *velx = 0; 
       *vely = 0; 
       *velz = 0;
       *press = pmin; 
       *eps = epsmin; 
       /* Bvec is already set above */
       /* What to do with eps???    */ 

       // Return here, circumventing poly_fallback completely //
       return fullstatus;

    }
 
    /* Choose when to allow a polytype EoS fallback */
    CCTK_INT poly_fallback = ( udens <= atmo_threshold ) ? 1 : 0; 

    /*************************************
     ** Call point-wise recovery scheme **
     *************************************/
    if ( !polytype ) {
        if ( utau > 1e100 ) {
           CCTK_VInfo(CCTK_THORNSTRING,"Poison in utau (%g), tau_in=%g, tau_min=%g  at (%g,%g,%g), level %d, weight %g, iteration %d",
                      utau, tau_in, tau_min, x, y, z, whisky_reflevel, carpetweight, iteration); 
        }
        if ( solve_system_method == HD_ONLY ) {
           ierr = EOS_solveHDsys( handle, false, rho_min, pmin,
                          udens, utau, usx, usy, usz,
                          gxx, gxy, gxz, gyy, gyz, gzz,
                          ugxx, ugxy, ugxz, ugyy, ugyz, ugzz,
                          &rho_new, &eps_new, &ux_new, &uy_new, &uz_new, soln_quality);
        } else {
           ierr = EOS_solvev2Z( handle, rho_min, atmo_tolerance, 
                          udens, utau, usx, usy, usz, *Bvecx, *Bvecy, *Bvecz,
                          gxx, gxy, gxz, gyy, gyz, gzz,
                          ugxx, ugxy, ugxz, ugyy, ugyz, ugzz,
                          &rho_new, &eps_new, &ux_new, &uy_new, &uz_new, soln_quality);
        }

        /*if ( ierr && B2 > MINB2 && (solve_system_method != MHD_5D_UiEps)) { // Try 5D
    		ierr = EOS_solvesys( handle, gstype, false, 
			udens, utau, usx, usy, usz, *Bvecx, *Bvecy, *Bvecz,
			gxx, gxy, gxz, gyy, gyz, gzz,
			ugxx, ugxy, ugxz, ugyy, ugyz, ugzz,
			&ux_new, &uy_new, &uz_new, &eps_new);
        }*/

	if ( ierr ) {
                fullstatus++;
		generalEOS_unsat=1;
                if ( poly_fallback || soln_quality[2] > 0 ) {poly_fallback = 1;}
		if ( (verbose > 2 && carpetweight > 0 && soln_quality[0]>0 ) || verbose > 3 ) 
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"GeneralEOS_solvesys returned an error at "
				"(x,y,z)=(%g,%g,%g) with weight %g. Errs=(%g,%g,%g). %s polytype.",
				x, y, z, carpetweight, soln_quality[0], soln_quality[1], soln_quality[2],
                                poly_fallback>0 ? "Trying" : "Not trying" );
		if ( verbose > 3 ) {
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     metric: g_ij=(%g,%g,%g,%g,%g,%g)", 
				gxx, gxy, gxz, gyy, gyz, gzz);
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     Cons: dens=%g, s_i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",dens_in,sx_in,sy_in,
				sz_in,Bnx_in,Bny_in,Bnz_in,tau_in);
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     Prims: rho=%g, eps=%g, v^i=(%g,%g,%g), B^i=(%g,%g,%g)",rho_in,eps_in,velx_in,
				vely_in,velz_in,Bvx_in,Bvy_in,Bvz_in);
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     Returned values: u^i = (%g,%g,%g), eps_new=%g",ux_new,uy_new,uz_new,eps_new);
		}
	}
	if ( eps_new < 0 ) {
                fullstatus++;
		generalEOS_unsat=1;
                if ( poly_fallback || soln_quality[2] > 0 ) {poly_fallback = 1;}
		if ( (verbose > 2 && carpetweight > 0 && soln_quality[0]>0) || verbose > 3  ) {
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"GeneralEOS_solvesys returned a negative internal energy at"
				"(x,y,z)=(%g,%g,%g) with weight %g.  Err=(%g,%g,%g). %s polytype.",
				x, y, z, carpetweight, soln_quality[0], soln_quality[1], soln_quality[2],
                                poly_fallback>0 ? "Trying" : "Not trying" );
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     metric: g_ij=(%g,%g,%g,%g,%g,%g)", 
				gxx, gxy, gxz, gyy, gyz, gzz);
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     Cons: dens=%g, s_i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",dens_in,sx_in,sy_in,
				sz_in,Bnx_in,Bny_in,Bnz_in,tau_in);
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     Prims: rho=%g, eps=%g, v^i=(%g,%g,%g), B^i=(%g,%g,%g)",rho_in,eps_in,velx_in,
				vely_in,velz_in,Bvx_in,Bvy_in,Bvz_in);
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				"     Returned values: u^i = (%g,%g,%g), eps_new=%g",ux_new,uy_new,uz_new,eps_new);
                }
	}

	/* Attempt to recover rest of primitives */
	if ( !generalEOS_unsat ) {

		*w_lorentz = sqrt( 1. + SQR(ux_new)*ugxx + SQR(uy_new)*ugyy + SQR(uz_new)*ugzz
			+ 2.*ux_new*uy_new*ugxy + 2.*ux_new*uz_new*ugxz + 2.*uy_new*uz_new*ugyz );
		*velx = ( ugxx*ux_new + ugxy*uy_new + ugxz*uz_new ) / *w_lorentz;
		*vely = ( ugxy*ux_new + ugyy*uy_new + ugyz*uz_new ) / *w_lorentz;
		*velz = ( ugxz*ux_new + ugyz*uy_new + ugzz*uz_new ) / *w_lorentz;

		//*rho = udens / *w_lorentz;
		*rho = rho_new; 
		if ( *rho < atmo_threshold ) {
			*rho = rho_min;
                        soln_quality[2] = soln_quality[2]+1.;
                        if ( carpetweight>=0.5 && !store_recovery_fails )
                           CCTK_VWarn(1,__LINE__,__FILE__, CCTK_THORNSTRING,
                                      "PHYS2ATMO! Set (%g,%g,%g), weight=%g to atmo in RecoverMHD after successful recovery yields rho<rho_min.",x,y,z,carpetweight);
			/* Handle the MHD atmosphere just like HD, leaving Bfield untouched. */
		}

		*eps = eps_new;
		*press = EOS_Pressure(handle, *rho, *eps );

                poly_fallback = 0;
	}

        // No polytype fallbacks for udens > atmo_threshold now that we're modifying the conservatives to deal with possible trouble spots!
        if ( generalEOS_unsat && carpetweight == 1 && soln_quality[0] > 0 && poly_fallback==0 ) {
           CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                      "WARNING! No solution found for r=(%g,%g,%g) on itn %d, level %d, with soln quality (%g,%g,%g) (eps_new=%g) from solve method %d. Cons_changed=%d",
                      x,y,z, iteration, whisky_reflevel, soln_quality[0], soln_quality[1], soln_quality[2], eps_new, solve_system_method, cons_have_changed);
           if ( verbose > 1 ) { /* Continue warning */
              /* End Plasma beta parameter */
              CCTK_REAL wL = sqrt( 1. + SQR(ux_new)*ugxx + SQR(uy_new)*ugyy + SQR(uz_new)*ugzz
			+ 2.*ux_new*uy_new*ugxy + 2.*ux_new*uz_new*ugxz + 2.*uy_new*uz_new*ugyz );
              rho_new = udens/wL;
              CCTK_REAL b2_comoving_out = ( B2 + SQR( (ux_new)*(*Bvecx) + (uy_new)*(*Bvecy) + (uz_new)*(*Bvecz) ) )/SQR(wL);
              CCTK_REAL betaPlasma_out = 0.5*b2_comoving_out/EOS_Pressure(handle, rho_new, eps_new);
              CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "         Incoming prims had beta_p = %15.7g (%15.7g), b2=%15.7g (%15.7g), B2=%15.7g, rho=%15.7g (%15.7g)",
                         betaPlasma_in, betaPlasma_out, b2_comoving_in, b2_comoving_out, B2, *rho, rho_new); 
           }
           
        }

    } // End general EOS attempts. Only do a polytype fallback if in atmo 

    //if ( polytype || udens < atmo_threshold ) {
    if ( polytype || poly_fallback ) { // Polytype atmosphere

        if ( solve_system_method == HD_ONLY ) {
           ierr = EOS_solveHDsys( phandle, true, rho_min, pmin,
                          udens, utau, usx, usy, usz,
                          gxx, gxy, gxz, gyy, gyz, gzz,
                          ugxx, ugxy, ugxz, ugyy, ugyz, ugzz,
                          &rho_new, &eps_new, &ux_new, &uy_new, &uz_new, soln_quality);

        } else {
           ierr = EOS_solveW( phandle, rho_min, atmo_tolerance,
                          udens, usx, usy, usz, *Bvecx, *Bvecy, *Bvecz,
                          gxx, gxy, gxz, gyy, gyz, gzz,
                          ugxx, ugxy, ugxz, ugyy, ugyz, ugzz,
                          &rho_new, &eps_new, &ux_new, &uy_new, &uz_new, soln_quality);

        }

        /*if ( ierr && B2 > MINB2 && (solve_system_method != MHD_5D_UiEps)) { // Try general (5D) method as last ditch
    		ierr = EOS_solvesys( phandle, gstype, true, 
			udens, utau, usx, usy, usz, *Bvecx, *Bvecy, *Bvecz,
			gxx, gxy, gxz, gyy, gyz, gzz,
			ugxx, ugxy, ugxz, ugyy, ugyz, ugzz,
			&ux_new, &uy_new, &uz_new, NULL, soln_quality);
        } */

	if ( ierr ) {
                if ( polytype ) fullstatus++; // We haven't registered a failure yet. Don't register fallback polytype attempts.
                if ( carpetweight < 1 && verbose > 2 ) {
			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                            "Not aborting code due to error in polytropeEOS_solvesys because point unphysical at (%g,%g,%g).",x,y,z);
                }
                if ( carpetweight == 1 ) {
		        if ( verbose > 2 ) {
			        CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				        "polytropeEOS_solvesys returned an error at (x,y,z)=(%g,%g,%g)"
				        "with weight %g. Err=(%g,%g,%g). Failed to find primitives given the following input:",
        				x, y, z, carpetweight, soln_quality[0], soln_quality[1], soln_quality[2] );
	        		CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
		        		"     metric: g_ij=(%g,%g,%g,%g,%g,%g)",
			        	gxx, gxy, gxz, gyy, gyz, gzz);
        			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
	        			"     Cons: dens=%g, s_i=(%g,%g,%g), B_n^i=(%g,%g,%g)",dens_in,sx_in,sy_in,
		        		sz_in,Bnx_in,Bny_in,Bnz_in);
			        if ( !polytype ) CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, "            tau=%g",tau_in);
        			CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
	        			"     Prims: rho=%g, eps=%g, v^i=(%g,%g,%g), B^i=(%g,%g,%g)",rho_in,eps_in,velx_in,
		        		vely_in,velz_in,Bvx_in,Bvy_in,Bvz_in);
			        CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
				        "     Returned values: u^i = (%g,%g,%g)",ux_new,uy_new,uz_new);
		        }
			//CCTK_WARN(0,"Aborting code due to error in polytropeEOS_solvesys on physical point: No solution for primitives found.");
			CCTK_WARN(1,"This run should fail shortly: No solution found for primitives even with polytropeEOS_solvesys.");
		}
	}

	*w_lorentz = sqrt( 1. + SQR(ux_new)*ugxx + SQR(uy_new)*ugyy + SQR(uz_new)*ugzz
		+ 2.*ux_new*uy_new*ugxy + 2.*ux_new*uz_new*ugxz + 2.*uy_new*uz_new*ugyz );
	*velx = ( ugxx*ux_new + ugxy*uy_new + ugxz*uz_new ) / *w_lorentz;
	*vely = ( ugxy*ux_new + ugyy*uy_new + ugyz*uz_new ) / *w_lorentz;
	*velz = ( ugxz*ux_new + ugyz*uy_new + ugzz*uz_new ) / *w_lorentz;

	*rho = udens / *w_lorentz;
	if ( *rho < atmo_threshold ) {
		*rho = rho_min;
                soln_quality[2] = soln_quality[2] + 1.;
                if ( carpetweight>=0.5 && !store_recovery_fails )
                   CCTK_VWarn(1,__LINE__,__FILE__, CCTK_THORNSTRING,
                      "PHYS2ATMO! Set (%g,%g,%g), weight=%g to atmo in RecoverMHD after polytrope recovery yields rho<rho_min.",x,y,z,carpetweight);
		/* Treat MHD atmosphere just like HD, not touching the Bfield */
	}

	*press = EOS_Pressure( phandle, *rho, JUNK );
        if ( eps != NULL ) {
           *eps = EOS_SpecificIntEnergy( phandle, *rho, JUNK );
        }

    } // END polytype attempt

    /* Extra physicality checks if point is physical, ignoring points where we've reset to atmosphere. Nothing done about these points here anyways. */
    if ( carpetweight==1 && soln_quality[2]==0 ) {

       /* Warn if W unphysical for physical points */
       if ( std::isnan(*w_lorentz) || *w_lorentz < 1) {
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"Con2PrimMHD: Unphysical Lorentz factor %g on carpet level %d, itn %d, r=(%g,%g,%g) after Con2Prim returns quality (%g,%g,%g)",
   		*w_lorentz, whisky_reflevel,iteration,x,y,z,soln_quality[0], soln_quality[1], soln_quality[2]);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     metric: g_ij=(%g,%g,%g,%g,%g,%g)",
   		gxx, gxy, gxz, gyy, gyz, gzz);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     Cons: dens=%g, s_i=(%g,%g,%g), B_n^i=(%g,%g,%g)",*dens,*sx,*sy,*sz,*Bnx,*Bny,*Bnz);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, "            tau=%g",*tau);
             CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                   "     Prims: rho=%g, v^i=(%g,%g,%g), B^i=(%g,%g,%g)",*rho,*velx,*vely,*velz,*Bvecx,*Bvecy,*Bvecz);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, "            eps=%g",*eps);
          CCTK_WARN(0,"Aborting code due to unphysical W for a physical point.");
       }
   
       /* Warn about speed if above warn_if_speed_above */
       CCTK_REAL speed2 = SQR(*velx)*gxx + SQR(*vely)*gyy + SQR(*velz)*gzz
                   + 2.*(*velx)*(*vely)*gxy + 2.*(*velx)*(*velz)*gxz + 2.*(*vely)*(*velz)*gyz;
       if ( speed2 > SQR(warn_if_speed_above) ) {
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"Con2PrimMHD: Speed of %g is above %g on carpet level %d after Con2Prim returns quality (%g,%g,%g)",
   		sqrt(speed2), warn_if_speed_above, whisky_reflevel, soln_quality[0], soln_quality[1], soln_quality[2]);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     Velocities = (%g,%g,%g)", *velx, *vely, *velz);	
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     Iteration %d at location (%g,%g,%g)", iteration, x, y, z);
       }
   
   
       /* Warn if eps is still negative and point is physical */
       if ( *eps < 0 && whisky_reflevel >= warn_from_reflevel ) {
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     Negative specific internal energy on carpet level %d, itn %d, r=(%g,%g,%g) after Con2Prim returns quality (%g,%g,%g)",
   		whisky_reflevel,iteration,x,y,z,soln_quality[0], soln_quality[1], soln_quality[2]);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     metric: g_ij=(%g,%g,%g,%g,%g,%g)",
   		gxx, gxy, gxz, gyy, gyz, gzz);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     Cons: dens=%g, s_i=(%g,%g,%g), B_n^i=(%g,%g,%g)",*dens,*sx,*sy,*sz,*Bnx,*Bny,*Bnz);
          CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, "            tau=%g",*tau);
   	  CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
   		"     Prims: rho=%g, v^i=(%g,%g,%g), B^i=(%g,%g,%g)",*rho,*velx,*vely,*velz,*Bvecx,*Bvecy,*Bvecz);
   	  CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING, "            eps=%g",*eps);
   	  CCTK_WARN(0,"Aborting code due to negative eps for a physical point.");
       }
    } else if ( *eps < 0 && verbose > 2 ) {
         CCTK_VWarn(2,__LINE__,__FILE__,CCTK_THORNSTRING,
		"     Negative specific internal energy, but I was told to ignore carpet level %d",whisky_reflevel);
    }

    return fullstatus;
}
