/* Solver function for the {Z=rho*h*W, v2} method of Noble et al 2006 */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "EOS_Base.h"
#include <stdio.h>

/* Local solver function */
static void Get_v2Z_fdf(const CCTK_REAL x[], struct solver_consts_hv2 *p, 
  CCTK_REAL dx[], CCTK_REAL *f, CCTK_REAL *df, 
  CCTK_REAL Jac[][2], CCTK_REAL resids[]);

/* Define parameter structure */
struct solver_consts_hv2 { 
    CCTK_INT handle, poly_switch;
    CCTK_REAL udens, utau, gamma_eos;
    CCTK_REAL uS2, uBS2, uB2; // S_i S^i, (B^i S_i)^2, and B_i B^i
    CCTK_REAL tau_plus_dens;
};

/***************************/
/* Primary solver routines */
/***************************/

CCTK_INT EOS_solvev2Z (CCTK_INT eos_handle, CCTK_REAL rho_min, CCTK_REAL atmo_tol,
    CCTK_REAL udens, CCTK_REAL utau, CCTK_REAL usx, CCTK_REAL usy, CCTK_REAL usz,
    CCTK_REAL bnx, CCTK_REAL bny, CCTK_REAL bnz,
    CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
    CCTK_REAL guxx, CCTK_REAL guxy, CCTK_REAL guxz, CCTK_REAL guyy, CCTK_REAL guyz, CCTK_REAL guzz,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *ux, CCTK_REAL *uy, CCTK_REAL *uz, CCTK_REAL soln_quality[3])
{

   DECLARE_CCTK_PARAMETERS;

   if ( verbose > 6 ) {
        CCTK_VInfo(CCTK_THORNSTRING,"Entered generalEOS_solvesys with the following quantities ..." );
        CCTK_VInfo(CCTK_THORNSTRING,"    g_ij = (%g,%g,%g,%g,%g,%g)",gxx,gxy,gxz,gyy,gyz,gzz);
        CCTK_VInfo(CCTK_THORNSTRING,"    g^ij = (%g,%g,%g,%g,%g,%g)",guxx,guxy,guxz,guyy,guyz,guzz);
   }

   /* Set parameters */ 
   struct solver_consts_hv2 param;
   CCTK_REAL uS2, uBS, uBS2, uB2, tau_plus_dens, eos_gamma;
   uS2 = SQR(usx)*guxx + SQR(usy)*guyy + SQR(usz)*guzz
         + 2.*usx*usy*guxy + 2.*usx*usz*guxz + 2.*usy*usz*guyz;
   uBS = bnx * usx + bny * usy + bnz * usz;
   uBS2 = SQR(uBS);
   uB2 = SQR(bnx)*gxx + SQR(bny)*gyy + SQR(bnz)*gzz
         + 2.*bnx*bny*gxy + 2.*bnx*bnz*gxz + 2.*bny*bnz*gyz;
   tau_plus_dens = utau + udens; 
   eos_gamma = EOS_Pressure(eos_handle, 1.0, 1.0) + 1.0; // HACK -- HARDCODED GAMMA LAW

   param.handle = eos_handle;
   param.poly_switch = 0;
   param.udens = udens;
   param.utau = utau;
   param.gamma_eos = eos_gamma;
   param.uS2 = uS2;
   param.uBS2 = uBS2;
   param.uB2 = uB2;
   param.tau_plus_dens = tau_plus_dens;

   /* used in reconstructing u_i */
   const CCTK_REAL Blowx = gxx*bnx + gxy*bny + gxz*bnz;
   const CCTK_REAL Blowy = gxy*bnx + gyy*bny + gyz*bnz;
   const CCTK_REAL Blowz = gxz*bnx + gyz*bny + gzz*bnz;

   /* Initial values */
   CCTK_REAL w_lorentz = sqrt( 1. + SQR(*ux)*guxx + SQR(*uy)*guyy + SQR(*uz)*guzz
                        + 2.**ux**uy*guxy + 2.**ux**uz*guxz + 2.**uy**uz*guyz );
   CCTK_REAL v2 = 1.-1./SQR(w_lorentz);
   CCTK_REAL trho = udens/w_lorentz, teps = *eps;
   CCTK_REAL tpress = EOS_Pressure(eos_handle, trho, teps);
   CCTK_REAL Z = (trho*(1+teps) + tpress)*SQR(w_lorentz);

   /* Check for poison */
   CCTK_INT poisoned[8];
   poisoned[0] = ( udens > 1e100 ) ? 1 : 0;
   poisoned[1] = ( utau > 1e100 ) ? 1 : 0;
   poisoned[2] = ( usx > 1e100 ) ? 1 : 0;
   poisoned[3] = ( usy > 1e100 ) ? 1 : 0;
   poisoned[4] = ( usz > 1e100 ) ? 1 : 0;
   poisoned[5] = ( bnx > 1e100 ) ? 1 : 0;
   poisoned[6] = ( bny > 1e100 ) ? 1 : 0;
   poisoned[7] = ( bnz > 1e100 ) ? 1 : 0;
   CCTK_INT anypoison = ( udens>1e100 || utau>1e100 || usx>1e100 || usy>1e100 
                          || usz>1e100 || bnx>1e100 || bny>1e100 || bnz>1e100 ) ? 1 : 0;

   /* Find Safe initial values (Cerda-Duran et al 2008) */
   CCTK_REAL epsmax = (utau - 0.5*uB2)/udens;
   CCTK_REAL rhomax = udens;
   CCTK_REAL Wmax = 100.;
   CCTK_REAL pressmax = EOS_Pressure( eos_handle, rhomax, epsmax );
   CCTK_REAL Zmax = utau + pressmax + udens - 0.5*uB2;
   CCTK_REAL Zmin = 0; //udens;
   CCTK_REAL v2max = 1. - 1./SQR(Wmax);

   /* Iterator */
   CCTK_INT itn = 0;
   CCTK_INT status = 1;
   CCTK_INT reset = 0;
   CCTK_INT fullcounter=0;
   CCTK_INT initial_success_at_itn=-1;
   CCTK_REAL local_tol=solver_tolerance;

   /* Reset info (DEBUG) */
   CCTK_INT Z_too_large = 0;
   CCTK_INT Z_too_small = 0;
   CCTK_INT v2_too_large = 0;
   CCTK_INT v2_too_small = 0;
   CCTK_INT hit_countmax = 0; 

   /* Set initial state */
   CCTK_REAL x_init[2] = { Z,v2 };
   CCTK_REAL x_old[2]  = { Z,v2 };
   CCTK_REAL x[2]      = { Z,v2 };
   CCTK_REAL dx[2];
   CCTK_REAL Jac[2][2];
   CCTK_REAL resids[2];
   CCTK_REAL f = 1.;
   CCTK_REAL df = 1.;
   CCTK_REAL errx = 1.;

   if ( verbose > 6 ) {
      CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",
                                       udens,usx,usy,usz,bnx,bny,bnz,utau);
      CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: v2 = %g Z=%g",x_init[0],x_init[1]);
   }

    
   CCTK_INT keep_going = 1;
   while (keep_going) {

        /* Get iteration info  */
        itn++;
        fullcounter++;
        Get_v2Z_fdf(x, &param, dx,&f,&df, Jac, resids);

        /* Store old */
        x_old[0] = x[0]; 
        x_old[1] = x[1]; 

        /* Take a step */
        x[0] += dx[0];
        x[1] += dx[1];

        /* Convergence criterion on Z */
        errx = (x[0]==0.) ? fabs(dx[0]) : fabs( dx[0]/x[0]);
        
        /* Impose physicality (with buffer) */
        if ( x[0] < Zmin ) { 
           Z_too_small++;
           x[0] = 0.75*Zmax;
           errx = 1.;
        } else {
           if ( x[0] > 2*Zmax ) {
              Z_too_large++;
              //CCTK_VInfo(CCTK_THORNSTRING," Z too large: %g -> %g ( > %g)", x_old[0], x[0], Zmax);
              x[0] = 0.75*Zmax;
              errx = 1.;
           }
        }
        if ( x[1] < 0. ) {
           x[1] = fabs(x[1]);
           v2_too_small++;
           errx = 1.;
        } else {
           if (x[1] > v2max) {
              x[1] = x_old[1];
              v2_too_large++;
              errx = 1.;
           }
        }

        /* Checks and Loop control */

        /* We've found a solution of enough accuracy */
        if ( fabs(errx) <= local_tol && (initial_success_at_itn<0) ) {
           initial_success_at_itn = itn;
           status = 0;
           keep_going = 0;
        }

        /* Switch keep_going back on if polishing the root */
        if ( status==0 ) {
           if ( (extra_solver_iterations>0) && itn<(initial_success_at_itn+extra_solver_iterations) ) {
               keep_going = 1;
           } else {
               keep_going = 0;
           }
        }

        /* Handle failure */
        if ( status>0 && (itn>=countmax) ) {
           hit_countmax++;
           /* 3rd failure, just quit trying */
           if ( reset > 1 ) {
              keep_going = 0;
           }
           /* 2nd failure, increase tolerance */
           if ( reset==1 ) {
              reset++;
              itn=0;
              x[0] = Zmax;
              x[1] = v2max;
              local_tol *= 100.;
           }
           /* 1st failure, set to Font's "safe" values from Cerda-Duran et al. 2008 */
           if ( reset==0 ) {
              reset++;
              itn=0;
              x[0] = Zmax;
              x[1] = v2max;
           }
        }

   }

   /* Get final measure of convergence */
   soln_quality[0] = errx;
   soln_quality[1] = fabs(f);

   /* Put new values into arguments only if we succeeded */
   if ( status == 0 ) {

      Z  = x[0]; 
      v2 = x[1];

      CCTK_REAL invW = sqrt(1.-v2);
      w_lorentz = 1./invW;
      trho = udens*invW;
      tpress =( (eos_gamma-1.)/eos_gamma )*( Z*SQR(invW) - trho ); // P(z,rho,gamma)
      //teps = ( Z*SQR(invW) - (tpress+trho) )/trho;
      teps = EOS_SpecificIntEnergy(eos_handle, trho, tpress); 
      *eps = teps;

      //if ( verbose > 2 )
      //   CCTK_VInfo(CCTK_THORNSTRING,"\tP,eps check: P=%g, eps=%g, Pfromeps=%g",
      //              tpress, teps, epscheck);

      if ( trho <= rho_min*(1.+atmo_tol) ) {
         *rho = rho_min;
         *ux = 0.;
         *uy = 0.;
         *uz = 0.;
         soln_quality[2] += 1;
      } else {
         *rho = trho;
         CCTK_REAL ufact = w_lorentz/(Z*(uB2+Z));
         *ux  = ufact * ( usx*Z + uBS*Blowx );
         *uy  = ufact * ( usy*Z + uBS*Blowy );
         *uz  = ufact * ( usz*Z + uBS*Blowz );
         //*ux=(Blowx*uBS+usx*Z)/(Z*(uB2+Z))*w_lorentz;
         //*uy=(Blowy*uBS+usy*Z)/(Z*(uB2+Z))*w_lorentz;
         //*uz=(Blowz*uBS+usz*Z)/(Z*(uB2+Z))*w_lorentz;
      }

       
      if ( verbose > 6 ) 
         CCTK_VInfo(CCTK_THORNSTRING,"\tSuccessful recovery: (v2,Z)=(%g,%g), u_i=(%g,%g,%g), eps=%g",
                    v2, (eps==NULL?-1:Z),*ux, *uy, *uz, (eps==NULL?-1:*eps) );
   }
  

   //if (verbose > 4 ) 
   if (verbose > 3 || ( verbose > 2 && status != 0 ) ) {
        char poisonline[35];
        sprintf(poisonline,"PRE-POISONED! (%d,%d,%d,%d,%d,%d,%d,%d)", poisoned[0], poisoned[1], poisoned[2], poisoned[3],
                poisoned[4], poisoned[5], poisoned[6], poisoned[7]);
        CCTK_VInfo(CCTK_THORNSTRING,"generalEOS solver: status = %d (%s) after %d (%d) itns. Reset=%d (%d-%d-%d-%d-%d), soln_quality=(%g,%g,%g), local_tol=%g. %s",
                   status, status ? "unsuccessful" : "successful" , itn, fullcounter, reset, hit_countmax, Z_too_small, Z_too_large,
                   v2_too_small, v2_too_large, soln_quality[0], soln_quality[1], soln_quality[2], local_tol, anypoison > 0 ? poisonline : " "  );
   }

   return status;

}

/********************************************/
/* Function to provide f, df, and residuals */
/********************************************/

static void Get_v2Z_fdf(const CCTK_REAL x[], struct solver_consts_hv2 *p, 
  CCTK_REAL dx[], CCTK_REAL *f, CCTK_REAL *df, 
  CCTK_REAL Jac[][2], CCTK_REAL resids[])
{

   /* Store current state */
   CCTK_REAL Z, v2;
   Z = x[0];
   v2 = x[1];

   /* Locals of parameters for now */
   CCTK_REAL dens, S2, BS2, B2, egamma, tau_plus_dens;
   dens = p->udens;
   //tau = p->utau;
   S2 = p->uS2;
   BS2 = p->uBS2;
   B2 = p->uB2;
   tau_plus_dens = p->tau_plus_dens;
   egamma = p->gamma_eos;

   /* Create useful temps, calculate derivs */
   CCTK_REAL invZ, W, invW, BdotSbyZ2, B2pZ, sqr_B2pZ;
   invW = sqrt(1.-v2);
   W = 1./invW;
   invZ = 1./Z;
   BdotSbyZ2 = BS2*SQR(invZ);
   B2pZ = Z+B2;
   sqr_B2pZ = SQR( B2pZ );

   CCTK_REAL gamfact, press, dpdv2, dpdZ;
   gamfact = (egamma-1.)/egamma;
   //rho = dens*invW;
   press = gamfact * (Z*SQR(invW)-dens*invW);
   dpdv2 = gamfact * (0.5*dens*W-Z);
   dpdZ  = gamfact * SQR(invW);

   /*****************************************************
    ** Function Evaluation                             **
    ** fS2 = S2 - v2*(B2+Z)^2 + (S.B)^2/Z^2*(B2+2Z)    **
    ** ftau = -Z - B2 + D + tau -(1/2)B2*( 1+v2 )      **
    **          + (1/2)(S.B)^2/Z^2                     **
    ** df = -( fS2^2 + ftau^2 )                        **
    ** f = 0.5*( df^2 )                                **
    *****************************************************/
   const CCTK_REAL fS2 = S2 - v2*sqr_B2pZ + (Z+B2pZ)*BdotSbyZ2;
   //const CCTK_REAL ftau = - Z - tau_plus_dens + press - 0.5*B2*invW + 0.5*BdotSbyZ2;
   const CCTK_REAL ftau = tau_plus_dens - Z + press + 0.5*BdotSbyZ2 - 0.5*B2*(1.+v2);

   resids[0] = fS2;
   resids[1] = ftau;
   *df = - SQR(fS2) - SQR(ftau);
   *f  = 0.5*SQR(*df);

   /*if ( std::isnan(*f) ) {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"df is NaN in v2Z solve: (fS2,ftau)=(%g,%g) "
         "(Z,v2)=(%g,%g), press=%g, tau+D=%g, S2=%g, B2=%g",
         fS2, ftau, Z, v2, press, tau_plus_dens, S2, B2 );
   }*/

   /*************************** 
    ** Set the Jacobian      ** 
    ** J[0][0] = d(fS2)/dZ   **
    ** J[0][1] = d(fS2)/dv2  **
    ** J[1][0] = d(ftau)/dZ  **
    ** J[1][1] = d(ftau)/dv2 **
    ***************************/

   CCTK_REAL dZ_fS2, dZ_ftau, dv2_fS2, dv2_ftau;
   dZ_fS2   = -2.*( B2pZ*v2 + BdotSbyZ2*( B2*invZ+1.) );
   dv2_fS2  = - sqr_B2pZ;
   dZ_ftau  = dpdZ - 1. - BdotSbyZ2*invZ;
   dv2_ftau = dpdv2 - 0.5*B2;
   Jac[0][0] = dZ_fS2; 
   Jac[0][1] = dv2_fS2;
   Jac[1][0] = dZ_ftau;
   Jac[1][1] = dv2_ftau;

  /***************************
   ** Calculate stepsize    **
   **   dx = (Jac^-1).f     **
   ***************************/

  //CCTK_REAL invdetJac = 1./( Jac[0][0]*Jac[1][1] - Jac[1][0]*Jac[0][1] );
  //dx[0] = - (   Jac[1][1]*fS2 - Jac[0][1]*ftau ) * invdetJac;
  //dx[1] = - ( - Jac[1][0]*fS2 + Jac[0][0]*ftau ) * invdetJac;
  //CCTK_REAL invdetJac = 1./( dZ_fS2*dv2_ftau - dv2_fS2*dZ_ftau );
  CCTK_REAL invdetJac = 1./( B2pZ*( B2pZ*dZ_ftau - 2.*dv2_ftau*(Z*v2 + BdotSbyZ2)*invZ ) ); 
  dx[0] = -(  dv2_ftau*fS2 - dv2_fS2*ftau ) * invdetJac;
  dx[1] = -( - dZ_ftau*fS2 +  dZ_fS2*ftau ) * invdetJac;

}
