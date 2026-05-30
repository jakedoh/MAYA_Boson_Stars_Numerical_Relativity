/* GSL-based solver functions for (u_i, epsilon) at a given point */
/* http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html */

#include "assert.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "EOS_Base.h"
#include "RecoverMHD.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <stdio.h>

#define CUBE(x) ((x)*(x)*(x))
#define EOS_DEpsByDRho(handle,rho,press) (-EOS_DPressByDRho(handle,rho,EOS_SpecificIntEnergy(handle,rho,press))/EOS_DPressByDEps(handle,rho,EOS_SpecificIntEnergy(handle,rho,press)))
#define EOS_DEpsByDPress(handle,rho,press) (1.0/EOS_DPressByDEps(handle,rho,EOS_SpecificIntEnergy(handle,rho,press)))

/* Local solver functions */
static int MHDsys_getv2Z_f(const gsl_vector * x, void * p, gsl_vector * f);
static int MHDsys_getv2Z_fdf(const gsl_vector * x, void * p, gsl_vector * f, gsl_matrix * J);
static int MHDsys_getv2Z_df(const gsl_vector * x, void * p, gsl_matrix * J);
static void print_state ( size_t iter, gsl_multiroot_fdfsolver * s);
static void print_state ( size_t iter, gsl_multiroot_fsolver * s);
static void get_final_errs ( const gsl_vector * dx, const gsl_vector * x, const int solver_n, double * max_rel_err, double * max_abs_err );

/* Parameter structure */
struct solver_params { 
    int handle, poly_switch;
    double udens, utau;
    double uS2, uBS2, uB2; // S_i S^i, (B^i S_i)^2, and B_i B^i
};

/***************************/
/* Primary solver routines */
/***************************/

int EOS_solvehv2sys (int eos_handle, enum GSL_SOLVER_TYPE stype,
    bool is_polytype,
    double udens, double utau, double usx, double usy, double usz,
    double bnx, double bny, double bnz,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double guxx, double guxy, double guxz, double guyy, double guyz, double guzz,
    double *ux, double *uy, double *uz, double *eps, double soln_quality[3])
{

   DECLARE_CCTK_PARAMETERS;

   if ( verbose > 6 ) {
        CCTK_VInfo(CCTK_THORNSTRING,"Entered %sEOS_solvesys with the following quantities ...", is_polytype ? "polytrope" : "general");
        CCTK_VInfo(CCTK_THORNSTRING,"    g_ij = (%g,%g,%g,%g,%g,%g)",gxx,gxy,gxz,gyy,gyz,gzz);
        CCTK_VInfo(CCTK_THORNSTRING,"    g^ij = (%g,%g,%g,%g,%g,%g)",guxx,guxy,guxz,guyy,guyz,guzz);
   }

   /* Choose solver type */
   bool is_deriv_solver=0;
   if ( stype < GST_HybridS )
        is_deriv_solver = 1;
   
   const size_t n = is_polytype ? 1 : 2;
   double uS2 = SQR(usx)*guxx + SQR(usy)*guyy + SQR(usz)*guzz
                + 2.*usx*usy*guxy + 2.*usx*usz*guxz + 2.*usy*usz*guyz;
   double uBS = bnx * usx + bny * usy + bnz * usz;
   double uBS2 = SQR(uBS);
   double uB2 = SQR(bnx)*gxx + SQR(bny)*gyy + SQR(bnz)*gzz
                + 2.*bnx*bny*gxy + 2.*bnx*bnz*gxz + 2.*bny*bnz*gyz;
   struct solver_params param = { eos_handle, is_polytype, 
         udens, is_polytype ? -1 : utau,
         uS2, uBS2, uB2,
         };

   /* used in reconstructing u_i */
   const double Blowx = gxx*bnx + gxy*bny + gxz*bnz;
   const double Blowy = gxy*bnx + gyy*bny + gyz*bnz;
   const double Blowz = gxz*bnx + gyz*bny + gzz*bnz;

   /* Solver type info */
   const gsl_multiroot_fsolver_type *Stype_f;
   const gsl_multiroot_fdfsolver_type *Stype_fdf;

   /* Hard-code Newton */
   Stype_f = NULL;
   Stype_fdf = gsl_multiroot_fdfsolver_newton;
   /*switch (stype)
   {
        case GST_Newton:
                Stype_f = NULL;
                Stype_fdf = gsl_multiroot_fdfsolver_newton;
                break;
        case GST_gNewton:
                Stype_f = NULL;
                Stype_fdf = gsl_multiroot_fdfsolver_gnewton;
                break;
        case GST_HybridSJ:
                Stype_f = NULL;
                Stype_fdf = gsl_multiroot_fdfsolver_hybridsj;
                break;
        case GST_HybridJ:
                Stype_f = NULL;
                Stype_fdf = gsl_multiroot_fdfsolver_hybridj;
                break;
        case GST_HybridS:
                Stype_f = gsl_multiroot_fsolver_hybrids;
                Stype_fdf = NULL;
                break;
        case GST_Hybrid:
                Stype_f = gsl_multiroot_fsolver_hybrid;
                Stype_fdf = NULL;
                break;
        case GST_dNewton:
                Stype_f = gsl_multiroot_fsolver_dnewton;
                Stype_fdf = NULL;
                break;
        case GST_Broyden:
                Stype_f = gsl_multiroot_fsolver_broyden;
                Stype_fdf = NULL;
                break;
        default:
                Stype_f = NULL;
                Stype_fdf = NULL;
                CCTK_WARN(0,"Invalid GSL Multiroot solver handle.");
                break;

   } */

   /* Set initial state */
   double w_lorentz = sqrt( 1. + SQR(*ux)*guxx + SQR(*uy)*guyy + SQR(*uz)*guzz
                        + 2.**ux**uy*guxy + 2.**ux**uz*guxz + 2.**uy**uz*guyz );
   double v2 = 1.-1./SQR(w_lorentz);
   double trho = udens/w_lorentz, teps = is_polytype ? EOS_SpecificIntEnergy(eos_handle, trho, 1.0) : *eps;
   double tpress = EOS_Pressure(eos_handle, trho, teps);
   double Z = (trho*(1+teps) + tpress)*SQR(w_lorentz);
   double x_init[2] = { v2, Z };
   gsl_vector_view x = gsl_vector_view_array(x_init, n);
   //printf("(udens, utau, S2, BS2, B2) = (%e %e %e %e %e)\n", udens, utau, uS2, uBS2, uB2); 

   /* Find Safe initial values (Cerda-Duran et al 2008) */
   double epsmax = (utau - 0.5*uB2)/udens;
   double rhomax = udens;
   double Wmax = 10000;
   double pressmax = EOS_Pressure( eos_handle, rhomax, epsmax );
   double Zmax = utau + pressmax + udens; // - 0.5*uB2;
   double Zmin = -udens; // udens
   double v2max = 1. - 1./SQR(Wmax);

   /* Check that initial guess is physical */
   if ( v2 > 1. || Z < 0 ) {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Unphysical initial guess ( v2, rho*h*W^2 )=(%g,%g)",v2,Z);
   }


   if ( verbose > 6 ) {
        if ( is_polytype ) {
                CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g)",
                                             udens,usx,usy,usz,bnx,bny,bnz);
                CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: v2 = %g",x_init[0]);
        } else {
                CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",
                                             udens,usx,usy,usz,bnx,bny,bnz,utau);
                CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: v2 = %g Z=%g",x_init[0],x_init[1]);
        }
   }

   /* Turn off internal gsl error handling so singular Jacobians 
      inside the puncture don't kill the simulation until more 
      useful information can be written. */
   gsl_set_error_handler_off();

   size_t iter = 0;
   int status = 0;
   int reset = 0;
   int loopcounter=0;
   double local_abstol=solver_epsabs;
   double local_reltol=solver_epsrel;
   if ( is_deriv_solver ) {
        /* Do entire iteration for deriv solvers */
        gsl_multiroot_function_fdf f = {MHDsys_getv2Z_f,
                                      MHDsys_getv2Z_df,
                                      MHDsys_getv2Z_fdf, n, &param};
        gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(Stype_fdf, n);
        gsl_multiroot_fdfsolver_set ( s, &f, &x.vector ); 

        if ( verbose > 5 )
                print_state (iter, s);

        int initial_success_at_iter = -1;
        do {
                iter++;
                gsl_vector_memcpy(&x.vector, s->x);
                status = gsl_multiroot_fdfsolver_iterate(s);

                /* If f or df contains an Inf or NaN (e.g. singular Jacobian) abort this. */
                if ( status == GSL_EBADFUNC )
                   break;                    

                /* Reset to Font's "safe" values if we get out of range */
                CCTK_REAL v2cur = gsl_vector_get(s->x,0);
                if ( n > 1 ) {
                   CCTK_REAL Zcur = gsl_vector_get (s->x,1);
                   if ( Zcur < Zmin || Zcur > Zmax || v2cur < 0. || v2cur > v2max ) { 
                       // Out of range. Reset to Font's "safe" values
                       if ( reset == 2 ) {
                          // We've already reset twice. Infinite loop? Kill it.
                          loopcounter += iter;
                          iter = countmax+extra_solver_iterations;
                       }
                       if ( reset == 1 ) {
                          // We've already reset once. Raise tolerance.
                          reset++;
                          loopcounter += iter;
                          local_abstol = 100.*local_abstol;
                          local_reltol = 100.*local_reltol;
                          gsl_vector_set(s->x, 0, v2max);
                          gsl_vector_set(s->x, 1, 0.75*Zmax);
                          gsl_vector_memcpy(&x.vector, s->x);
                          gsl_multiroot_fdfsolver_set ( s, &f, &x.vector );
                          iter = 0;
                       }
                       if ( reset == 0 ) {
                          reset++;
                          loopcounter += iter;
                          iter = 0;
                          gsl_vector_set(s->x, 0, v2max);
                          gsl_vector_set(s->x, 1, 0.75*Zmax);
                          gsl_vector_memcpy(&x.vector, s->x);
                          gsl_multiroot_fdfsolver_set ( s, &f, &x.vector );
                          //if ( verbose > 5 ) CCTK_WARN(1,"Resetting to Font since out of range.");
                       }
                   }
                } else {
                   if ( v2cur < 0. || v2cur > v2max ) { // v2 out of range. Reset to Font's "safe" value
                       // Out of range. Reset to Font's "safe" values
                       if ( reset==2 ) {
                          // We've already reset twice. Infinite loop? Kill it.
                          loopcounter += iter;
                          iter = countmax+extra_solver_iterations;
                       }
                       if ( reset==1 ) {
                          // We've already reset once. Raise tolerance.
                          reset++;
                          loopcounter += iter;
                          iter = 0;
                          local_abstol = 100.*local_abstol;
                          local_reltol = 100.*local_reltol;
                          gsl_vector_set(s->x, 0, v2max);
                          gsl_vector_memcpy(&x.vector, s->x);
                          gsl_multiroot_fdfsolver_set ( s, &f, &x.vector );
                       }
                       if ( reset==0 ) {
                          reset++;
                          loopcounter += iter;
                          iter = 0;
                          gsl_vector_set(s->x, 0, v2max);
                          gsl_vector_memcpy(&x.vector, s->x);
                          gsl_multiroot_fdfsolver_set ( s, &f, &x.vector );
                          //if ( verbose > 5 ) CCTK_WARN(1,"Resetting to Font since out of range.");
                       }
                   }
                }
                //if ( verbose > 5 )
                //        print_state (iter, s);

                if (status) /* Is the solver stuck or failed? */
                        break; 

                /* The convergence errors for this method depends most strongly on |dZ/Z|, 
                   so use only test_delta.  solver_epsabs should be very small. */
                status = gsl_multiroot_test_delta ( s->dx, s->x, local_abstol, local_reltol );

                /* Override GSL_SUCCESS for extra_solver_iterations iterations if we are not going to hit countmax */
                if ( status == GSL_SUCCESS ) {
                   if ( initial_success_at_iter == -1 )
                     initial_success_at_iter = iter;
                   if ( iter < initial_success_at_iter+extra_solver_iterations && iter < (size_t)countmax )
                     status = GSL_CONTINUE;
                }

                /* If we're about to hit countmax and we haven't reset, reset. */
                if ( reset<2 && iter == (size_t)countmax ) {
                   reset++;
                   if ( reset == 1 ) {
                      local_abstol=100.*local_abstol;
                      local_reltol=100.*local_reltol;
                   }
                   loopcounter += iter;
                   iter = 0;
                   gsl_vector_set(s->x, 0, v2max);
                   gsl_vector_set(s->x, 1, Zmax);
                   gsl_vector_memcpy(&x.vector, s->x);
                   gsl_multiroot_fdfsolver_set ( s, &f, &x.vector );
                }

        } while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

        /* Put new values into arguments only if we succeeded */
        if ( status == GSL_SUCCESS ) {
                double eos_gamma = EOS_Pressure(eos_handle, 1.0, 1.0) + 1.0; // HACK
                v2 = gsl_vector_get( s->x, 0);
                w_lorentz = 1/sqrt(1.-v2);
                trho = udens/w_lorentz;
                if ( !is_polytype ) {
                        Z = gsl_vector_get( s->x, 1);
                        // Assumes gamma-law EOS
                        teps = (Z-trho*SQR(w_lorentz))/(trho*SQR(w_lorentz)*eos_gamma);
                        *eps = teps;
                }
                // Equ. 47 of Giacomazzo
                *ux=(Blowx*uBS+usx*Z)/(Z*(uB2+Z))*w_lorentz;
                *uy=(Blowy*uBS+usy*Z)/(Z*(uB2+Z))*w_lorentz;
                *uz=(Blowz*uBS+usz*Z)/(Z*(uB2+Z))*w_lorentz;

                if ( verbose > 6 ) 
                   CCTK_VInfo(CCTK_THORNSTRING,"\tSuccessful recovery: (v2,Z)=(%g,%g), u_i=(%g,%g,%g), eps=%g",
                             v2, (eps==NULL?-1:Z),*ux, *uy, *uz, (eps==NULL?-1:*eps) );

        }

        /* Get last measure of convergence, regardless if successful */
        CCTK_REAL max_relerr = 0;
        CCTK_REAL max_abserr = 0;
        for ( int i=0; i<n; i++ ) {
           double dxL = fabs( gsl_vector_get( s->dx, i ) );
           double xL  = fabs( gsl_vector_get( s->x, i ) );
           double Vabserr = fabs( dxL - local_reltol*xL );
           double Vrelerr = fabs( dxL - local_abstol ) / xL;
           if ( Vabserr > max_abserr ) 
              max_abserr = Vabserr;
           if ( Vrelerr > max_relerr ) 
              max_relerr = Vrelerr;
        } 
        soln_quality[0] = max_relerr;
        soln_quality[1] = max_abserr;

        gsl_multiroot_fdfsolver_free(s);

   } else {
       CCTK_WARN(0,"GeneralEOS system_type \"h v2\" currently only supports derivative solvers. Please alter parameter solver_type.");
   }

   //if (verbose > 4 ) {
   if (verbose > 2 || ( verbose > 1 && status != GSL_SUCCESS ) ) {
        CCTK_VInfo(CCTK_THORNSTRING,"%sEOS solver: status = %d (%s) after %d (%d) iterations. Reset=%d, soln_quality=(%g,%g,%g)",is_polytype ? "polytrope" : "general",
                   status,gsl_strerror(status),(int)iter,loopcounter, reset, soln_quality[0], soln_quality[1], soln_quality[2] );
   }

   return status != GSL_SUCCESS;

}
/*******************************************/
/* Functions to print current solver state */
/*******************************************/

static void print_state ( size_t iter, gsl_multiroot_fdfsolver * s)
{
        DECLARE_CCTK_PARAMETERS;
        CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (v2,Z)=(% e % e) "
                "f(v2,Z) = (% e % e) J(v2,Z) = ((% e,% e),(% e,% e))", (int)iter,
                gsl_vector_get( s->x, 0 ), s->x->size < 2 ? -1 : gsl_vector_get( s->x, 1),
                gsl_vector_get( s->f, 0 ), s->f->size < 2 ? -1 : gsl_vector_get( s->f, 1),
                gsl_matrix_get( s->J, 0,0 ), s->f->size < 2 ? -1 : gsl_matrix_get( s->J, 0,1),
                s->f->size < 2 ? -1 : gsl_matrix_get( s->J, 1,0), 
                s->f->size < 2 ? -1 : gsl_matrix_get( s->J, 1,1)
                );
}

static void print_state ( size_t iter, gsl_multiroot_fsolver * s)
{
        CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (v2,Z)=(% .3f % .3f) "
                "f(v2,Z) = (% .3e % .3e)", (int)iter,
                gsl_vector_get( s->x, 0 ), s->x->size < 2 ? -1 : gsl_vector_get( s->x, 1),
                gsl_vector_get( s->f, 0 ), s->f->size < 2 ? -1 : gsl_vector_get( s->f, 1));
}

/***************************************/
/* Functions to provide f, df, and fdf */
/***************************************/

static int MHDsys_getv2Z_f(const gsl_vector * x, void * p, gsl_vector * f)
{
  return MHDsys_getv2Z_fdf(x, p, f, NULL);
}

static int MHDsys_getv2Z_fdf(const gsl_vector * x, void * p, gsl_vector * f, gsl_matrix * J)
{

   struct solver_params * param = (struct solver_params *)p;

   /* Get information from the parameters into local variables */
   int handle, do_polytrope;
   handle = (param->handle);
   do_polytrope = (param->poly_switch);

   double dens, tau, S2, BS2, B2;
   dens = (param->udens);
   tau = (param->utau);
   S2 = (param->uS2);
   BS2 = (param->uBS2);
   B2 = (param->uB2);

   /* Get current W, Z from vector x */
   /* Get pressure, its derivatives, and the specific enthalpy */
   double v2, Z, invZ, W, invW, rho, eps, press;
   v2 = gsl_vector_get(x,0);
   W = 1./sqrt(1-v2);
   invW = 1./W;

   rho = dens*invW;
   double eos_gamma = EOS_Pressure(handle, 1.0, 1.0) + 1.0; // HACK
   if ( do_polytrope == 0 ) {
        Z = gsl_vector_get(x,1);
        press = (eos_gamma-1)/eos_gamma * (Z*(1-v2)-rho);
   } else {
        press = EOS_Pressure(handle,rho,1.0);
        eps = EOS_SpecificIntEnergy(handle,rho,1.0);
        Z = (rho*(1+eps)+press)*SQR(W);
   }

   /* helper variables */
   invZ = 1./Z;

   /************************** 
    ** Function Evaluations ** 
    **************************/

   if ( NULL != f ) {
      /* Construct s_i definitions */
      const double fS2 = -1*(SQR(Z+B2) * v2 - (2*Z+B2)*BS2*SQR(invZ) - S2);
      gsl_vector_set (f, 0, fS2);

      /* Construct tau definition if not polytrope */
      if ( do_polytrope == 0 ) {
           const double ftau  = Z + B2 - press - 0.5*B2*(1-v2) - 0.5*BS2*SQR(invZ) - dens - tau;
           gsl_vector_set (f, 1, ftau);
      }
   }

   /************************** 
    ** Jacobian Evaluations ** 
    **************************/

   if ( NULL != J ) {
      /* Shorthands to make eqns readable */
      double dfS2dv2 = -SQR(Z+B2);
      gsl_matrix_set (J, 0, 0, dfS2dv2);

      /* press derivatives and derivatives of tau if not polytrope */
      if ( do_polytrope == 0 ) {
           //double dfS2dZ = 2*(Z+B2)*v2 - 2*BS2*SQR(invZ) + 2*(2*Z+B2)*BS2*CUBE(invZ);
           double dfS2dZ = -1*(2*(Z+B2)*v2 + 2.*BS2*( B2*CUBE(invZ) + SQR(invZ) ));
           gsl_matrix_set (J, 0, 1, dfS2dZ);

           double dpdv2 = (eos_gamma-1)/eos_gamma*( 0.5*dens*W - Z);
           double dpdZ = (eos_gamma-1)/eos_gamma*(1-v2);

           double dftaudv2 = 0.5*B2 - dpdv2;
           double dftaudZ = 1.0 + BS2*CUBE(invZ) - dpdZ;
           gsl_matrix_set (J, 1, 0, dftaudv2);
           gsl_matrix_set (J, 1, 1, dftaudZ);
      }
   }

   return GSL_SUCCESS;

}

static int MHDsys_getv2Z_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
   return MHDsys_getv2Z_fdf(x, p, NULL, J);
}
