/* GSL-based solver functions for (u_i, epsilon) at a given point */
/* http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html */
/* uses method of arXiv:1001.0575v2 (1d method works only for Gamma-law) */

#include "assert.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "EOS_Base.h"
#include "RecoverMHD.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <stdio.h>

#define CUBE(x) ((x)*(x)*(x))
#define EOS_DEpsByDRho(handle,rho,press) (-EOS_DPressByDRho(handle,rho,EOS_SpecificIntEnergy(handle,rho,press))/EOS_DPressByDEps(handle,rho,EOS_SpecificIntEnergy(handle,rho,press)))
#define EOS_DEpsByDPress(handle,rho,press) (1.0/EOS_DPressByDEps(handle,rho,EOS_SpecificIntEnergy(handle,rho,press)))

/* Local solver functions */
static double MHDsys_getUE_f(const double x, void * p);
static void MHDsys_getUE_fdf(const double x, void * p, double * y, double * dy);
static double MHDsys_getUE_df(const double x, void * p);
static void print_state ( size_t iter, gsl_root_fdfsolver * s);
static void print_state ( size_t iter, gsl_root_fsolver * s);
static void get_final_errs ( const double * dx, const double * x, const int solver_n, double * max_rel_err, double * max_abs_err );

/* Parameter structure */
struct solver_params { 
    int handle, poly_switch;
    double udens, utau;
    double uS2, uBS2, uB2; // S_i S^i, (B^i S_i)^2, and B_i B^i
};

/***************************/
/* Primary solver routines */
/***************************/

int EOS_solvehW2sys (int eos_handle, enum GSL_SOLVER_TYPE stype,
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
   
   const size_t n = 1;
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
   const gsl_root_fsolver_type *Stype_f;
   const gsl_root_fdfsolver_type *Stype_fdf;

   /* hack we misuse the names introduce for the multiroot solvers */
   switch (stype)
   {
	case GST_Newton:
		Stype_f = NULL;
		Stype_fdf = gsl_root_fdfsolver_newton;
		break;
	case GST_gNewton:
		Stype_f = NULL;
		Stype_fdf = gsl_root_fdfsolver_secant;
		break;
	case GST_HybridSJ:
		Stype_f = NULL;
		Stype_fdf = gsl_root_fdfsolver_steffenson;
		break;
	case GST_HybridS:
		Stype_f = gsl_root_fsolver_bisection;
		Stype_fdf = NULL;
		break;
	case GST_Hybrid:
		Stype_f = gsl_root_fsolver_falsepos;
		Stype_fdf = NULL;
		break;
	case GST_dNewton:
		Stype_f = gsl_root_fsolver_brent;
		Stype_fdf = NULL;
		break;
	case GST_HybridJ:
	case GST_Broyden:
	default:
		Stype_f = NULL;
		Stype_fdf = NULL;
		CCTK_WARN(0,"Invalid GSL Multiroot solver handle.");
		break;

   }

   /* Set initial state */
   double w_lorentz = sqrt( 1. + SQR(*ux)*guxx + SQR(*uy)*guyy + SQR(*uz)*guzz
			+ 2.**ux**uy*guxy + 2.**ux**uz*guxz + 2.**uy**uz*guyz );
   double rho = udens/w_lorentz;
   double eps_temp = is_polytype ? EOS_SpecificIntEnergy( eos_handle, rho, 1.0 ) : *eps;
   double press_temp =  EOS_Pressure( eos_handle, rho, eps_temp );
   double rho_enthalpy = rho * (1.+eps_temp) + press_temp;
   double UE_init[1] = { rho_enthalpy * SQR(w_lorentz) };
   gsl_vector_view x = gsl_vector_view_array(UE_init, n);

   // make sure the initial guess is physical
   double X0 = UE_init[0];
   while((SQR(X0)*SQR(X0 + uB2))/(SQR(X0)*SQR(X0+uB2)-SQR(X0)*uS2-(2*X0+uB2)*uBS2) < 1.) {
       X0 *= 10.;
   }
   UE_init[0] = X0;

   if ( verbose > 6 ) {
        if ( is_polytype ) {
		CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g)",
					     udens,usx,usy,usz,bnx,bny,bnz);
		CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: x = %g",UE_init[0]);
	} else {
		CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",
					     udens,usx,usy,usz,bnx,bny,bnz,utau);
		CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: x = %g",UE_init[0]);
	}
   }

   size_t iter = 0;
   int status;
   if ( is_deriv_solver ) {
	/* Do entire iteration for deriv solvers */
   	gsl_function_fdf f = {MHDsys_getUE_f,
			      MHDsys_getUE_df,
			      MHDsys_getUE_fdf, &param};
   	gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(Stype_fdf);
	gsl_root_fdfsolver_set ( s, &f, UE_init[0] ); 

	if ( verbose > 5 )
		print_state (iter, s);

        double x_old, x_new, dx;
	do {
		iter++;
                x_old = gsl_root_fdfsolver_root(s);
		status = gsl_root_fdfsolver_iterate(s);
                x_new =  gsl_root_fdfsolver_root(s);

                /* validate current guess and correct of out of bounds */
                dx = x_new-x_old;
                int i = 0;
                while((SQR(x_new)*SQR(x_new + uB2))/(SQR(x_new)*SQR(x_new+uB2)-SQR(x_new)*uS2-(2*x_new+uB2)*uBS2) < 1.0 && i < 10)
                {
                    x_new -= 0.1*dx;
                }
                if(x_new < 0 || i > 0)
	            gsl_root_fdfsolver_set ( s, &f, fabs(x_new) ); 

		if ( verbose > 5 )
                {
                        double rho, press, teps, tux,tuy,tuz, h, X;
                        X = gsl_root_fdfsolver_root( s );
                        w_lorentz = sqrt((SQR(X)*SQR(X + uB2))/(SQR(X)*SQR(X+uB2)-SQR(X)*uS2-(2*X+uB2)*uBS2));
                        rho = udens/w_lorentz;
                        if ( !is_polytype ) {
                                double eos_gamma;
                                /* hack to get eos_gamma out of the equation of state */
                                eos_gamma = EOS_Pressure(eos_handle, 1.0, 1.0) + 1.0;

                                /* Equ. 36 */
                                press = (eos_gamma-1.0)/eos_gamma * (X/SQR(w_lorentz) - udens/w_lorentz);

                                teps=EOS_SpecificIntEnergy(eos_handle,rho,press);
                        }
                        // Equ. 47 of Giacomazzo (x = Z)
                        tux=(Blowx*uBS+usx*X)/(X*(uB2+X))*w_lorentz;
                        tuy=(Blowy*uBS+usy*X)/(X*(uB2+X))*w_lorentz;
                        tuz=(Blowz*uBS+usz*X)/(X*(uB2+X))*w_lorentz;
                        CCTK_VInfo(CCTK_THORNSTRING,"(ux uy uz eps)=(% e %e %e %e)",
                                tux,tuy,tuz,teps);

			print_state (iter, s);
                }

		if (status) /* Is the solver stuck or failed? */
			break; 

                if (test_residual)
                        status = gsl_root_test_residual ( GSL_FN_FDF_EVAL_F(&f, gsl_root_fdfsolver_root(s)), solver_tolerance );
                else
                        status = gsl_root_test_delta ( x_old, x_new, solver_epsabs, solver_epsrel );
	} while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

   	/* Put new values into arguments */
   	if ( status == GSL_SUCCESS ) {
                double rho, press, teps, h, X;
                X = gsl_root_fdfsolver_root( s );
                w_lorentz = sqrt((SQR(X)*SQR(X + uB2))/(SQR(X)*SQR(X+uB2)-SQR(X)*uS2-(2*X+uB2)*uBS2));
                rho = udens/w_lorentz;
		if ( !is_polytype ) {
                        double eos_gamma;
                        /* hack to get eos_gamma out of the equation of state */
                        eos_gamma = EOS_Pressure(eos_handle, 1.0, 1.0) + 1.0;

                        /* Equ. 36 */
                        press = (eos_gamma-1.0)/eos_gamma * (X/SQR(w_lorentz) - udens/w_lorentz);

			teps=EOS_SpecificIntEnergy(eos_handle,rho,press);
                        *eps = teps;
                }
                // Equ. 47 of Giacomazzo (x = Z)
		*ux=(Blowx*uBS+usx*X)/(X*(uB2+X))*w_lorentz;
		*uy=(Blowy*uBS+usy*X)/(X*(uB2+X))*w_lorentz;
		*uz=(Blowz*uBS+usz*X)/(X*(uB2+X))*w_lorentz;
	}

        /* Get last measure of convergence, regardless if successful */
        soln_quality[0] = ( fabs(dx) - solver_epsabs )/fabs(x_new); 
        soln_quality[1] = ( fabs(dx) - solver_epsrel*fabs(x_new) );

   	gsl_root_fdfsolver_free(s);

   } else {
	/* Do entire iteration for non-deriv solvers */
   	gsl_function f = {MHDsys_getUE_f, &param};
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(Stype_f);
	gsl_root_fsolver_set ( s, &f, 0.5*UE_init[0], 1.5*UE_init[0] ); 

	if ( verbose > 5 )
		print_state (iter, s);

        double x_old, x_new, dx;
	do {
		iter++;
                x_old = gsl_root_fsolver_root(s);
		status = gsl_root_fsolver_iterate(s);
                x_new =  gsl_root_fsolver_root(s);

                /* Store this */
                dx = x_new-x_old;

		if ( verbose > 5 )
			print_state (iter, s);

		if (status) /* Is the solver stuck or failed? */
			break; 

                if (test_residual)
                        status = gsl_root_test_residual ( GSL_FN_EVAL(&f, gsl_root_fsolver_root(s)), solver_tolerance );
                else
                        status = gsl_root_test_delta ( x_old, x_new, solver_epsabs, solver_epsrel );
	} while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

   	/* Put new values into arguments */
   	if ( status == GSL_SUCCESS ) {
                double rho, press, teps, h, X;
                X = gsl_root_fsolver_root( s );
                w_lorentz = sqrt((SQR(X)*SQR(X + uB2))/(SQR(X)*SQR(X+uB2)-SQR(X)*uS2-(2*X+uB2)*uBS2));
                rho = udens/w_lorentz;
		if ( !is_polytype ) {
                        double eos_gamma;
                        /* hack to get eos_gamma out of the equation of state */
                        eos_gamma = EOS_Pressure(eos_handle, 1.0, 1.0) + 1.0;

                        /* Equ. 36 */
                        press = (eos_gamma-1.0)/eos_gamma * (X/SQR(w_lorentz) - udens/w_lorentz);

			teps=EOS_SpecificIntEnergy(eos_handle,rho,press);
                        *eps = teps;
                }
                // Equ. 47 of Giacomazzo (x = Z)
		*ux=(Blowx*uBS+usx*X)/(X*(uB2+X))*w_lorentz;
		*uy=(Blowy*uBS+usy*X)/(X*(uB2+X))*w_lorentz;
		*uz=(Blowz*uBS+usz*X)/(X*(uB2+X))*w_lorentz;
	}

        /* Get last measure of convergence, regardless if successful */
        soln_quality[0] = ( fabs(dx) - solver_epsabs )/fabs(x_new); 
        soln_quality[1] = ( fabs(dx) - solver_epsrel*fabs(x_new) );

   	gsl_root_fsolver_free(s);

   }

   if (verbose > 4 ) {
	CCTK_VInfo(CCTK_THORNSTRING,"GeneralEOS solver: status = %d (%s) after %d iterations.",status,gsl_strerror(status),(int)iter);
   }

   return status != GSL_SUCCESS;

}
/*******************************************/
/* Functions to print current solver state */
/*******************************************/

static void print_state ( size_t iter, gsl_root_fdfsolver * s)
{
	CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (x)=% e",
		(int)iter,
		gsl_root_fdfsolver_root( s ));
}

static void print_state ( size_t iter, gsl_root_fsolver * s)
{
	CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (x)=% e",
		(int)iter,
		gsl_root_fsolver_root( s ));
}

/***************************************/
/* Functions to provide f, df, and fdf */
/***************************************/

static double MHDsys_getUE_f(const double x, void * p)
{
  double retval;
  MHDsys_getUE_fdf(x, p, &retval, NULL);
  return retval;
}

static void MHDsys_getUE_fdf(const double X, void * p, double * y, double * dy)
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
   double W2, W, invW, press, eos_gamma;
   W2 = (SQR(X)*SQR(X + B2))/(SQR(X)*SQR(X+B2)-SQR(X)*S2-(2*X+B2)*BS2);
   W = W2 > 1. ? sqrt(W2) : 1.0;
   invW = 1/W;

   /* hack to get eos_gamma out of the equation of state */
   eos_gamma = EOS_Pressure(handle, 1.0, 1.0) + 1.0;

   /* Equ. 36 */
   press = (eos_gamma-1.0)/eos_gamma * (X*SQR(invW) - dens*invW);

   /************************** 
    ** Function Evaluations ** 
    **************************/

   if ( NULL != y ) {
      /* Equ. 34 */
      const double f0 = X-press-0.5*(BS2/SQR(X)+B2*SQR(invW)) + B2 - tau - dens;
      *y = f0;
   }

   /************************** 
    ** Jacobian Evaluations ** 
    **************************/

   if ( NULL != dy ) {
      /* Shorthands to make eqns readable */
       /* all from Maple */
      double dW2dx =2*X*SQR(X+B2)/(SQR(X)*SQR(X+B2)-SQR(X)*S2-(2*X+B2)*BS2)+2*SQR(X)*(X+B2)/(SQR(X)*SQR(X+B2)-SQR(X)*S2-(2*X+B2)*BS2)-SQR(X)*SQR(X+B2)*(2*X*SQR(X+B2)+2*SQR(X)*(X+B2)-2*X*S2-2*BS2)/SQR(SQR(X)*SQR(X+B2)-SQR(X)*S2-(2*X+B2)*BS2); 
      double dpressdx = (eos_gamma-1.0)*SQR(invW)/eos_gamma + (eos_gamma-1.0)/eos_gamma * (-2*X*CUBE(invW) + dens*SQR(invW)) * dW2dx/(2*W);
      double df0dx = 1.0-dpressdx-0.5*(-3*BS2/CUBE(X)-B2*SQR(SQR(invW))*dW2dx);

      *dy = df0dx;
   }

   return;

}

static double MHDsys_getUE_df(const double x, void * p)
{
   double retval;
   MHDsys_getUE_fdf(x, p, NULL, &retval);
   return retval;
}
