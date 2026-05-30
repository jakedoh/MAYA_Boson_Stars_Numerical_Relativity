/* Solver function for the {Z=rho*h*W, v2} method of Noble et al 2006 */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "EOS_Base.h"
#include <stdio.h>

#define CUBE(x) ((x)*(x)*(x))

/* Local solver function */
static void Get_W_fdf(const CCTK_REAL x, struct solver_consts_W *p, 
    CCTK_REAL *dx, CCTK_REAL *f, CCTK_REAL *df ); 

/* Define parameter structure */
struct solver_consts_W { 
    CCTK_INT handle, poly_switch;
    CCTK_REAL udens, gamma_eos;
    CCTK_REAL uS2, uBS2, uB2; // S_i S^i, (B^i S_i)^2, and B_i B^i
};

/***************************/
/* Primary solver routines */
/***************************/

CCTK_INT EOS_solveW (CCTK_INT eos_handle, CCTK_REAL rho_min, CCTK_REAL atmo_tol,
    CCTK_REAL udens, CCTK_REAL usx, CCTK_REAL usy, CCTK_REAL usz,
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
   struct solver_consts_W param;
   CCTK_REAL uS2, uBS, uBS2, uB2, eos_gamma;
   uS2 = SQR(usx)*guxx + SQR(usy)*guyy + SQR(usz)*guzz
         + 2.*usx*usy*guxy + 2.*usx*usz*guxz + 2.*usy*usz*guyz;
   uBS = bnx * usx + bny * usy + bnz * usz;
   uBS2 = SQR(uBS);
   uB2 = SQR(bnx)*gxx + SQR(bny)*gyy + SQR(bnz)*gzz
         + 2.*bnx*bny*gxy + 2.*bnx*bnz*gxz + 2.*bny*bnz*gyz;
   eos_gamma = EOS_Pressure(eos_handle, 1.0, 1.0) + 1.0; // HACK -- HARDCODED GAMMA LAW

   param.handle = eos_handle;
   param.poly_switch = 1;
   param.udens = udens;
   param.gamma_eos = eos_gamma;
   param.uS2 = uS2;
   param.uBS2 = uBS2;
   param.uB2 = uB2;

   /* used in reconstructing u_i */
   const CCTK_REAL Blowx = gxx*bnx + gxy*bny + gxz*bnz;
   const CCTK_REAL Blowy = gxy*bnx + gyy*bny + gyz*bnz;
   const CCTK_REAL Blowz = gxz*bnx + gyz*bny + gzz*bnz;


   /* Find Safe initial values */
   CCTK_REAL Wmax = 100.;

   /* Initial values */
   CCTK_REAL w_lorentz_in = sqrt( 1. + SQR(*ux)*guxx + SQR(*uy)*guyy + SQR(*uz)*guzz
                              + 2.**ux**uy*guxy + 2.**ux**uz*guxz + 2.**uy**uz*guyz );
   CCTK_REAL w_lorentz = w_lorentz_in > Wmax ? Wmax : w_lorentz_in ;


   /* Iterator */
   CCTK_INT itn = 0;
   CCTK_INT status = 1;
   CCTK_INT reset = 0;
   CCTK_INT fullcounter=0;
   CCTK_INT initial_success_at_itn=-1;
   CCTK_REAL local_tol=solver_tolerance;

   /* Set initial state */
   CCTK_REAL x_init = w_lorentz_in;
   CCTK_REAL x_old  = w_lorentz_in;
   CCTK_REAL x = w_lorentz;
   CCTK_REAL dx = 1;;
   CCTK_REAL f = 1.;
   CCTK_REAL df = 1.;
   CCTK_REAL errx = 1.;

   if ( verbose > 6 ) {
      CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g)",
                                       udens,usx,usy,usz,bnx,bny,bnz);
      CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: W = %g",x_init);
   }

    
   CCTK_INT keep_going = 1;
   while (keep_going) {

        /* Get iteration info  */
        itn++;
        fullcounter++;
        Get_W_fdf(x, &param, &dx, &f, &df );

        /* Store old */
        x_old = x; 

        /* Take a step */
        x += dx;

        /* Convergence criterion on Z */
        CCTK_REAL Z, Z_old;
        Z     = udens*(1. + eos_gamma*EOS_SpecificIntEnergy(eos_handle, udens/x,1.0)) * x;
        Z_old = udens*(1. + eos_gamma*EOS_SpecificIntEnergy(eos_handle, udens/x_old,1.0)) * x_old;
        errx = Z==0. ? fabs( Z-Z_old) : fabs( (Z-Z_old)/Z );
 
        /* Is this physical? Start at Wmax if not */
        if ( x < 1. || x > Wmax ) { 
           x = Wmax;
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
        if ( status>0 && (itn>=countmax) && reset==0 ) {
           /* First failure, set to Font's "safe" values from Cerda-Duran et al. 2008 */
           if ( reset==0 ) {
              reset++;
              itn=0;
              x = Wmax;
           }
           /* 2nd failure, increase tolerance */
           if ( reset==1 ) {
              reset++;
              itn=0;
              x = Wmax;
              local_tol *= 100.;
           }
           /* 3rd failure, just quit trying */
           if ( reset > 1 ) {
              keep_going = 0;
           }
        }

   }

   /* Put new values into arguments only if we succeeded */
   if ( status == 0 ) {

      w_lorentz = x; 

      CCTK_REAL invW, enthalpy, temp_rho, temp_eps;
      invW = 1./w_lorentz;
      temp_rho = udens*invW;
      temp_eps = EOS_SpecificIntEnergy(eos_handle,temp_rho,1.0);
      //*eps = temp_eps;

      if ( temp_rho <= rho_min*(1.+atmo_tol) ) {
         *rho = rho_min;
         *ux = 0.;
         *uy = 0.;
         *uz = 0.;
         soln_quality[2] += 1;
      } else {
         *rho = temp_rho;
         CCTK_REAL Z, enthalpy, ufact;
         enthalpy = 1. + eos_gamma*temp_eps;
         Z = udens*enthalpy*w_lorentz;
         ufact = w_lorentz/(Z*(uB2+Z*w_lorentz));
         *ux  = ufact * ( usx*Z + uBS*Blowx );
         *uy  = ufact * ( usy*Z + uBS*Blowy );
         *uz  = ufact * ( usz*Z + uBS*Blowz );
      }

      //*ux=(Blowx*uBS+usx*Z)/(Z*(uB2+Z))*w_lorentz;
      //*uy=(Blowy*uBS+usy*Z)/(Z*(uB2+Z))*w_lorentz;
      //*uz=(Blowz*uBS+usz*Z)/(Z*(uB2+Z))*w_lorentz;
       
      if ( verbose > 6 ) 
         CCTK_VInfo(CCTK_THORNSTRING,"\tSuccessful recovery: W=%g, u_i=(%g,%g,%g), eps=%g",
                    w_lorentz, *ux, *uy, *uz, temp_eps );
   }
  

   /* Get final measure of convergence */
   soln_quality[0] = errx;
   soln_quality[1] = fabs(f);

   //if (verbose > 4 ) 
   if (verbose > 2 || ( verbose > 1 && status != 0 ) ) {
        CCTK_VInfo(CCTK_THORNSTRING,"polytropeEOS solver: status = %d (%s) after %d (%d) itns. Reset=%d, soln_quality=(%g,%g,%g)",
                   status, status ? "unsuccessful" : "successful" , itn, fullcounter, reset, soln_quality[0], soln_quality[1], soln_quality[2] );
   }

   return status;

}

/********************************************/
/* Function to provide f, df, and residuals */
/********************************************/

static void Get_W_fdf(const CCTK_REAL x, struct solver_consts_W *p, 
                      CCTK_REAL *dx, CCTK_REAL *f, CCTK_REAL *df ) 
{

   /* Store current state */
   CCTK_REAL W = x;

   /* Locals of parameters for now */
   CCTK_REAL dens, S2, BS2, B2, egamma;
   dens = p->udens;
   S2 = p->uS2;
   BS2 = p->uBS2;
   B2 = p->uB2;
   egamma = p->gamma_eos;
   CCTK_INT handle = p->handle;

   /* Create useful temps, calculate derivs */
   CCTK_REAL invW, rho, press, h, dhdrho;
   invW = 1./W;
   rho = dens*invW;
   press = EOS_Pressure(handle,rho,1.0);
   h = 1. + (egamma-1.)*EOS_SpecificIntEnergy(handle,rho,1.0);
   dhdrho = egamma*press/SQR(rho);

   CCTK_REAL Z, invZ, BdotSbyZ2, B2pZ, sqr_B2pZ;
   Z = dens*h*W;
   invZ = 1./Z;
   BdotSbyZ2 = BS2*SQR(invZ);
   B2pZ = Z+B2;
   sqr_B2pZ = SQR( B2pZ );

   /*****************************************************
    ** Function Evaluation                             **
    ** fS2 = S2 - v2*(B2+Z)^2 + (S.B)^2/Z^2*(B2+2Z)    **
    ** df = d(fS2)/dW                                  **
    *****************************************************/
   const CCTK_REAL fS2 = sqr_B2pZ * (1-SQR(invW)) - (Z + B2pZ)*BdotSbyZ2 - S2;
   double dZdW = dens*h - SQR(dens)*invW*dhdrho;
   double dS2dZ = 2*(Z+B2)*(1-SQR(invW)) + 2*BS2*SQR(invZ) + 2*B2*BS2*CUBE(invZ);
   double dW_fS2 = dS2dZ * dZdW + 2*CUBE(invW)*SQR(Z+B2);

   //CCTK_REAL dZdW = dens*h - SQR(dens)*invW*dhdrho;
   //CCTK_REAL dS2dZ = 2*B2pZ*(1-SQR(invW)) + 2*BdotSbyZ2 + 2*B2*BdotSbyZ2*invZ;
   //CCTK_REAL dW_fS2 = dS2dZ * dZdW + 2*CUBE(invW)*sqr_B2pZ;

   *df = dW_fS2;
   *f  = fS2;
   *dx = - (*f)/(*df);

   if ( std::isnan(*f) ) {
      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"NaN in fS2. Z=%g, B2pZ=%g, S2=%g, invW=%g, BdotSbyZ2=%g",
         Z, B2pZ, S2, invW, BdotSbyZ2 );
   }
   if ( std::isnan(x + *dx) ) {
      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"NaN in next step. W->W'=%g->%g, dx=%g, f=%g, df=%g, Z=%g, dZdW=%g, dS2dZ=%g, dhdrho=%g, dens=%g, rho=%g",
         W, W+(*dx), (*dx), fS2, dW_fS2, Z, dZdW, dS2dZ, dhdrho, dens, rho );
   }

}
