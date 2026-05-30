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
#include <gsl/gsl_multiroots.h>
#include <stdio.h>

#define CUBE(x) ((x)*(x)*(x))
#define EOS_DEpsByDRho(handle,rho,press) (-EOS_DPressByDRho(handle,rho,EOS_SpecificIntEnergy(handle,rho,press))/EOS_DPressByDEps(handle,rho,EOS_SpecificIntEnergy(handle,rho,press)))
#define EOS_DEpsByDPress(handle,rho,press) (1.0/EOS_DPressByDEps(handle,rho,EOS_SpecificIntEnergy(handle,rho,press)))

/* Local solver functions */
static int MHDsys_getUE_f(const gsl_vector * x, void * p, gsl_vector * f);
static int MHDsys_getUE_fdf(const gsl_vector * x, void * p, gsl_vector * f, gsl_matrix * J);
static int MHDsys_getUE_df(const gsl_vector * x, void * p, gsl_matrix * J);
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

int EOS_solveWPresssys (int eos_handle, enum GSL_SOLVER_TYPE stype,
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
   double uB2 = SQR(bnx)*guxx + SQR(bny)*guyy + SQR(bnz)*guzz
		+ 2.*bnx*bny*guxy + 2.*bnx*bnz*guxz + 2.*bny*bnz*guyz;
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

   switch (stype)
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

   }

   /* Set initial state */
   double w_lorentz = sqrt( 1. + SQR(*ux)*guxx + SQR(*uy)*guyy + SQR(*uz)*guzz
			+ 2.**ux**uy*guxy + 2.**ux**uz*guxz + 2.**uy**uz*guyz );
   double UE_init[2] = { w_lorentz, is_polytype ? -1 : EOS_Pressure(eos_handle, udens/w_lorentz, *eps) };
   gsl_vector_view x = gsl_vector_view_array(UE_init, n);
   printf("(udens, utau, S2, BS2, B2) = (%e %e %e %e %e)\n", udens, utau, uS2, uBS2, uB2); 

   if ( verbose > 6 ) {
        if ( is_polytype ) {
		CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g)",
					     udens,usx,usy,usz,bnx,bny,bnz);
		CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: W = %g",UE_init[0]);
	} else {
		CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",
					     udens,usx,usy,usz,bnx,bny,bnz,utau);
		CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: W = %g press=%g",UE_init[0],UE_init[1]);
	}
   }

   size_t iter = 0;
   int status;
   if ( is_deriv_solver ) {
	/* Do entire iteration for deriv solvers */
   	gsl_multiroot_function_fdf f = {MHDsys_getUE_f,
				      MHDsys_getUE_df,
				      MHDsys_getUE_fdf, n, &param};
   	gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(Stype_fdf, n);
	gsl_multiroot_fdfsolver_set ( s, &f, &x.vector ); 

	if ( verbose > 5 )
		print_state (iter, s);

	do {
		iter++;
                gsl_vector_memcpy(&x.vector, s->x);
		status = gsl_multiroot_fdfsolver_iterate(s);

                bool reset = false;
                if ( gsl_vector_get (s->x, 0) < 1.0) {
                    gsl_vector_set(s->x, 0, 1.0);
                    reset = true;
                }
#define W_TOO_BIG 1e9
                if ( gsl_vector_get( s->x, 0) > W_TOO_BIG ) {
                    gsl_vector_set(s->x, 0, gsl_vector_get (&x.vector, 0));
                    reset = true;
                }
                if ( gsl_vector_get (s->x, 1) < 0.0) {
                    gsl_vector_set(s->x, 1, 1e-15);
                    reset = true;
                }
                if(reset) {
                    gsl_vector_memcpy(&x.vector, s->x);
	            gsl_multiroot_fdfsolver_set ( s, &f, &x.vector ); 
                }

		if ( verbose > 5 )
			print_state (iter, s);

		if (status) /* Is the solver stuck or failed? */
			break; 

                if (test_residual)
                        status = gsl_multiroot_test_residual ( s->f, solver_tolerance );
                else
                        status = gsl_multiroot_test_delta ( s->dx, s->x, solver_epsabs, solver_epsrel );
	} while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

   	/* Put new values into arguments */
   	if ( status == GSL_SUCCESS ) {
                double rho, press, teps, h, Z;
                w_lorentz = gsl_vector_get( s->x, 0);
                rho = udens/w_lorentz;
		if ( !is_polytype ) {
                        press=gsl_vector_get( s->x, 1);
			teps=EOS_SpecificIntEnergy(eos_handle,rho,press);
                        *eps = teps;
                } else {
                        press=EOS_Pressure(eos_handle,rho,1.0);
			teps=EOS_SpecificIntEnergy(eos_handle,rho,1.0);
                }
                h = 1 + teps + press/rho;
                Z = udens*h*w_lorentz;
                // Equ. 47
		*ux=(Blowx*uBS+usx*Z)/(Z*(uB2+Z))*w_lorentz;
		*uy=(Blowy*uBS+usy*Z)/(Z*(uB2+Z))*w_lorentz;
		*uz=(Blowz*uBS+usz*Z)/(Z*(uB2+Z))*w_lorentz;
	}

        /* Get last measure of convergence, regardless if successful */
        CCTK_REAL max_relerr, max_abserr;
        get_final_errs ( s->dx, s->x, n, &max_relerr, &max_abserr );
        soln_quality[0] = max_relerr;
        soln_quality[1] = max_abserr;

   	gsl_multiroot_fdfsolver_free(s);

   } else {
	/* Do entire iteration for non-deriv solvers */
   	gsl_multiroot_function f = {MHDsys_getUE_f, n, &param};
	gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(Stype_f, n);
	gsl_multiroot_fsolver_set ( s, &f, &x.vector ); 

	if ( verbose > 5 )
		print_state (iter, s);

	do {
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		if ( verbose > 5 )
			print_state (iter, s);

		if (status != GSL_SUCCESS) { /* Is the solver stuck or failed? */
                        if (status == GSL_ENOPROG) {
                            status = GSL_CONTINUE;
                            continue;
                        }
                        if (status == GSL_ENOPROG)
                                status = gsl_multiroot_test_residual ( s->f, solver_tolerance );
			break; 
                }

                if (test_residual)
                        status = gsl_multiroot_test_residual ( s->f, solver_tolerance );
                else
                        status = gsl_multiroot_test_delta ( s->dx, s->x, solver_epsabs, solver_epsrel );
	} while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

   	/* Put new values into arguments */
   	if ( status == GSL_SUCCESS ) {
                double rho, press, teps, h, Z;
                w_lorentz = gsl_vector_get( s->x, 0);
                rho = udens/w_lorentz;
		if ( !is_polytype ) {
                        press=gsl_vector_get( s->x, 1);
			teps=EOS_SpecificIntEnergy(eos_handle,rho,press);
                        *eps = teps;
                } else {
                        press=EOS_SpecificIntEnergy(eos_handle,rho,1.0);
			teps=EOS_SpecificIntEnergy(eos_handle,rho,1.0);
                        *eps = teps;
                }
                h = 1 + teps + press/rho;
                Z = udens*h*w_lorentz;
                // Equ. 47
		*ux=(Blowx*uBS+usx*Z)/(Z*(uB2+Z))*w_lorentz;
		*uy=(Blowy*uBS+usy*Z)/(Z*(uB2+Z))*w_lorentz;
		*uz=(Blowz*uBS+usz*Z)/(Z*(uB2+Z))*w_lorentz;
	}

        /* Get last measure of convergence, regardless if successful */
        CCTK_REAL max_relerr, max_abserr;
        get_final_errs ( s->dx, s->x, n, &max_relerr, &max_abserr );
        soln_quality[0] = max_relerr;
        soln_quality[1] = max_abserr;

   	gsl_multiroot_fsolver_free(s);

   }

   if (verbose > 4 ) {
	CCTK_VInfo(CCTK_THORNSTRING,"GeneralEOS solver: status = %d (%s) after %d iterations.",status,gsl_strerror(status),(int)iter);
   }

   // Returns 0 if success. 1 Otherwise
   return status != GSL_SUCCESS;

}
/*******************************************/
/* Functions to print current solver state */
/*******************************************/

static void print_state ( size_t iter, gsl_multiroot_fdfsolver * s)
{
	CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (W,press)=(% e % e) "
		"f(W,press) = (% e % e) J(W,press) = ((% e,% e),(% e,% e))", (int)iter,
		gsl_vector_get( s->x, 0 ), s->x->size < 2 ? -1 : gsl_vector_get( s->x, 1),
		gsl_vector_get( s->f, 0 ), s->f->size < 2 ? -1 : gsl_vector_get( s->f, 1),
		gsl_matrix_get( s->J, 0,0 ), s->f->size < 2 ? -1 : gsl_matrix_get( s->J, 0,1),
                s->f->size < 2 ? -1 : gsl_matrix_get( s->J, 1,0), 
                s->f->size < 2 ? -1 : gsl_matrix_get( s->J, 1,1)
                );
}

static void print_state ( size_t iter, gsl_multiroot_fsolver * s)
{
	CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (W,press)=(% .3f % .3f) "
		"f(W,press) = (% .3e % .3e)", (int)iter,
		gsl_vector_get( s->x, 0 ), s->x->size < 2 ? -1 : gsl_vector_get( s->x, 1),
		gsl_vector_get( s->f, 0 ), s->f->size < 2 ? -1 : gsl_vector_get( s->f, 1));
}

static void get_final_errs ( const gsl_vector * dx, const gsl_vector * x, const int solver_n, double * max_rel_err, double * max_abs_err )
{

       DECLARE_CCTK_PARAMETERS;

       *max_rel_err=0;
       *max_abs_err=0;
       for ( int i=0; i<solver_n; i++ ) {
           double dxL = fabs( gsl_vector_get( dx, i ) );
           double xL  = fabs( gsl_vector_get(  x, i ) );
           double local_abserr = dxL - solver_epsrel*xL;
           double local_relerr = ( dxL - solver_epsabs ) / xL;
           if ( local_abserr > *max_abs_err ) 
              *max_abs_err = local_abserr;
           if ( local_relerr > *max_rel_err ) 
              *max_rel_err = local_relerr;
       } 

}

/***************************************/
/* Functions to provide f, df, and fdf */
/***************************************/

static int MHDsys_getUE_f(const gsl_vector * x, void * p, gsl_vector * f)
{
  return MHDsys_getUE_fdf(x, p, f, NULL);
}

static int MHDsys_getUE_fdf(const gsl_vector * x, void * p, gsl_vector * f, gsl_matrix * J)
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

   /* Get current W, press from vector x */
   /* Get pressure, its derivatives, and the specific enthalpy */
   double W, invW, rho, eps, press, dhdp, dhdrho, h;
   W = gsl_vector_get(x,0);
   invW = 1/W;

   rho = dens*invW;
   if ( do_polytrope == 0 ) {
        press = gsl_vector_get(x,1);
   } else {
        press = EOS_Pressure(handle,rho,1.0);
   }
   eps = EOS_SpecificIntEnergy(handle,rho,press);
   dhdp = EOS_DEpsByDPress(handle,rho,press) + 1/rho;
   dhdrho = EOS_DEpsByDRho(handle,rho,press) - press/SQR(rho);
   h = 1 + eps + press/rho;

   /* helper variables */
   double Z = dens * h * W;
   double invZ = 1/Z;

   /************************** 
    ** Function Evaluations ** 
    **************************/

   if ( NULL != f ) {
      /* Construct s_i definitions */
      const double fS2 = SQR(Z+B2) * (1-SQR(invW)) - (2*Z+B2)*BS2*SQR(invZ) - S2;
      gsl_vector_set (f, 0, fS2);

      /* Construct tau definition if not polytrope */
      if ( do_polytrope == 0 ) {
           const double ftau  = Z + B2 - press - B2/(2*SQR(W)) - BS2/(2*SQR(Z)) - dens - tau;
           gsl_vector_set (f, 1, ftau);
      }
   }

   /************************** 
    ** Jacobian Evaluations ** 
    **************************/

   if ( NULL != J ) {
      /* Shorthands to make eqns readable */
      double dZdW = dens*h - SQR(dens)*invW*dhdrho;
      double dS2dZ = 2*(Z+B2)*(1-SQR(invW)) + 2*BS2*SQR(invZ) + 2*B2*BS2*CUBE(invZ);
      double dfS2dW = dS2dZ * dZdW + 2*CUBE(invW)*SQR(Z+B2);
      double dfS2dp = dS2dZ * dhdp * dens*W;

      gsl_matrix_set (J, 0, 0, dfS2dW);

      /* press derivatives and derivatives of tau if not polytrope */
      if ( do_polytrope == 0 ) {
           gsl_matrix_set (J, 0, 1, dfS2dp);

           double dZdp = dens*W*dhdp;
           double dftaudW = (1+BS2*CUBE(invZ))*dZdW + B2*CUBE(invW);
           double dftaudp = (1+BS2*CUBE(invZ))*dZdp - 1.;
           gsl_matrix_set (J, 1, 0, dftaudW);
           gsl_matrix_set (J, 1, 1, dftaudp);
      }
   }

   return GSL_SUCCESS;

}

static int MHDsys_getUE_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
   return MHDsys_getUE_fdf(x, p, NULL, J);
}
