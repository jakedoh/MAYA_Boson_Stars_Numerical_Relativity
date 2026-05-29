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
    double usx, usy, usz;
    double bvecx, bvecy, bvecz;
    double gxx, gxy, gxz, gyy, gyz, gzz;
    double guxx, guxy, guxz, guyy, guyz, guzz;
};

/***************************/
/* Primary solver routines */
/***************************/

int EOS_solvesys (int eos_handle, enum GSL_SOLVER_TYPE stype,
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
        if ( is_polytype ) {
		CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g)",
					     udens,usx,usy,usz,bnx,bny,bnz);
		CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: u_i = (%g,%g,%g)",*ux,*uy,*uz);
	} else {
		CCTK_VInfo(CCTK_THORNSTRING,"    Undensitized conservatives: dens=%g, s^i=(%g,%g,%g), B_n^i=(%g,%g,%g), tau=%g",
					     udens,usx,usy,usz,bnx,bny,bnz,utau);
		CCTK_VInfo(CCTK_THORNSTRING,"    Initial guess: u_i = (%g,%g,%g) eps=%g",*ux,*uy,*uz,*eps);
	}
   }

   /* Choose solver type */
   bool is_deriv_solver=0;
   if ( stype < GST_HybridS )
	is_deriv_solver = 1;
   
   const size_t n = is_polytype ? 3 : 4;
   struct solver_params param = { eos_handle, is_polytype, 
	 udens, is_polytype ? -1 : utau,
	 usx, usy, usz,
	 bnx, bny, bnz,
	 gxx, gxy, gxz, gyy, gyz, gzz,
	 guxx, guxy, guxz, guyy, guyz, guzz };

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
   double UE_init[4] = { *ux, *uy, *uz, is_polytype ? -1 : *eps };
   gsl_vector *x = gsl_vector_alloc(n);
   gsl_vector_set ( x, 0, UE_init[0] );
   gsl_vector_set ( x, 1, UE_init[1] );
   gsl_vector_set ( x, 2, UE_init[2] );
   if ( !is_polytype )
   	gsl_vector_set ( x, 3, UE_init[3] );

   size_t iter = 0;
   int status;
   if ( is_deriv_solver ) {
	/* Do entire iteration for deriv solvers */
   	gsl_multiroot_function_fdf f = {MHDsys_getUE_f,
				      MHDsys_getUE_df,
				      MHDsys_getUE_fdf, n, &param};
   	gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(Stype_fdf, n);
	gsl_multiroot_fdfsolver_set ( s, &f, x ); 

	if ( verbose > 5 )
		print_state (iter, s);

	do {
		iter++;
		status = gsl_multiroot_fdfsolver_iterate(s);

		if ( verbose > 5 )
			print_state (iter, s);

		if (status != GSL_SUCCESS) { /* Is the solver stuck or failed? */
                        if (status == GSL_ENOPROG) /* this triggers (at least) when the initial guess was already perfect */
                                status = gsl_multiroot_test_delta ( s->dx, s->x, solver_epsabs, solver_epsrel );
			break; 
                }

                if (test_residual)
                        status = gsl_multiroot_test_residual ( s->f, solver_tolerance );
                else
                        status = gsl_multiroot_test_delta ( s->dx, s->x, solver_epsabs, solver_epsrel );
	} while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

   	/* Put new values into arguments */
   	if ( status == GSL_SUCCESS ) {
		*ux=gsl_vector_get( s->x, 0);
		*uy=gsl_vector_get( s->x, 1);
		*uz=gsl_vector_get( s->x, 2);
		if ( !is_polytype )
			*eps=gsl_vector_get( s->x, 3);
	}

   	gsl_multiroot_fdfsolver_free(s);
   	gsl_vector_free(x);

   } else {
	/* Do entire iteration for non-deriv solvers */
   	gsl_multiroot_function f = {MHDsys_getUE_f, n, &param};
	gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(Stype_f, n);
	gsl_multiroot_fsolver_set ( s, &f, x ); 

	if ( verbose > 5 )
		print_state (iter, s);

	do {
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		if ( verbose > 5 )
			print_state (iter, s);

		if (status != GSL_SUCCESS) { /* Is the solver stuck or failed? */
                        if (status == GSL_ENOPROG)
                                status = gsl_multiroot_test_delta ( s->dx, s->x, solver_epsabs, solver_epsrel );
			break; 
                }

                if (test_residual)
                        status = gsl_multiroot_test_residual ( s->f, solver_tolerance );
                else
                        status = gsl_multiroot_test_delta ( s->dx, s->x, solver_epsabs, solver_epsrel );
	} while ( status == GSL_CONTINUE && iter < (size_t)countmax ); 

   	/* Put new values into arguments */
   	if ( status == GSL_SUCCESS ) {
		*ux=gsl_vector_get( s->x, 0);
		*uy=gsl_vector_get( s->x, 1);
		*uz=gsl_vector_get( s->x, 2);
		if ( !is_polytype )
			*eps=gsl_vector_get( s->x, 3);
	}

        /* Get last measure of convergence, regardless if successful */
        CCTK_REAL max_relerr, max_abserr;
        get_final_errs ( s->dx, s->x, n, &max_relerr, &max_abserr );
        soln_quality[0] = max_relerr;
        soln_quality[1] = max_abserr;

        /* Free it all */
   	gsl_multiroot_fsolver_free(s);
   	gsl_vector_free(x);

   }

   if (verbose > 4 ) {
	CCTK_VInfo(CCTK_THORNSTRING,"GeneralEOS solver: status = %s after %d iterations.",gsl_strerror(status),(int)iter);
   }

   return status;

}
/*******************************************/
/* Functions to print current solver state */
/*******************************************/

static void print_state ( size_t iter, gsl_multiroot_fdfsolver * s)
{
	CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (u_i,eps)=(% .3f % .3f % .3f % .3f) "
		"f(u_i,eps) = (% .3e % .3e % .3e % .3e)", (int)iter,
		gsl_vector_get( s->x, 0 ), gsl_vector_get( s->x, 1),
		gsl_vector_get( s->x, 2 ), s->x->size < 4 ? -1 : gsl_vector_get( s->x, 3),
		gsl_vector_get( s->f, 0 ), gsl_vector_get( s->f, 1),
		gsl_vector_get( s->f, 2 ), s->f->size < 4 ? -1 : gsl_vector_get( s->f, 3));
}

static void print_state ( size_t iter, gsl_multiroot_fsolver * s)
{
	CCTK_VInfo(CCTK_THORNSTRING,"iter=%d (u_i,eps)=(% .3f % .3f % .3f % .3f) "
		"f(u_i,eps) = (% .3e % .3e % .3e % .3e)", (int)iter,
		gsl_vector_get( s->x, 0 ), gsl_vector_get( s->x, 1),
		gsl_vector_get( s->x, 2 ), s->x->size < 4 ? -1 : gsl_vector_get( s->x, 3),
		gsl_vector_get( s->f, 0 ), gsl_vector_get( s->f, 1),
		gsl_vector_get( s->f, 2 ), s->f->size < 4 ? -1 : gsl_vector_get( s->f, 3));
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

   double dens, tau, sx, sy, sz, Bnx, Bny, Bnz;
   dens = (param->udens);
   tau = (param->utau);
   sx = (param->usx); /* lower */
   sy = (param->usy);
   sz = (param->usz);
   Bnx = (param->bvecx); /* upper */
   Bny = (param->bvecy);
   Bnz = (param->bvecz);

   double gxx, gxy, gxz, gyy, gyz, gzz;
   gxx = (param->gxx);
   gxy = (param->gxy);
   gxz = (param->gxz);
   gyy = (param->gyy);
   gyz = (param->gyz);
   gzz = (param->gzz);

   double guxx, guxy, guxz, guyy, guyz, guzz;
   guxx = (param->guxx);
   guxy = (param->guxy);
   guxz = (param->guxz);
   guyy = (param->guyy);
   guyz = (param->guyz);
   guzz = (param->guzz);

   /* Get current u_i, eps from vector x */
   double ux, uy, uz, eps;
   ux  = gsl_vector_get(x,0);
   uy  = gsl_vector_get(x,1);
   uz  = gsl_vector_get(x,2);
   if ( do_polytrope == 0 ) {
   	eps = gsl_vector_get(x,3);
   } else {
   	eps = -1;
   }

   /* Construct quantity w == alpha u^0 */
   const double w2 = 1. + SQR(ux)*guxx + SQR(uy)*guyy + SQR(uz)*guzz
		 + 2.*ux*uy*guxy + 2.*ux*uz*guxz + 2.*uy*uz*guyz;
   const double w = sqrt(w2);
   const double invw2 = 1/w2, invw = 1/w;

   /* construct v^i = gamma^{ij} u_j/W */
   const double velx = (guxx*ux + guxy*uy + guxz*uz) * invw;
   const double vely = (guxy*ux + guyy*uy + guyz*uz) * invw;
   const double velz = (guxz*ux + guyz*uy + guzz*uz) * invw;

   /* Construct b_i, and b^2 */
   const double Bnlowx = gxx*Bnx + gxy*Bny + gxz*Bnz;
   const double Bnlowy = gxy*Bnx + gyy*Bny + gyz*Bnz;
   const double Bnlowz = gxz*Bnx + gyz*Bny + gzz*Bnz;
               
   const double udotB = Bnx*ux + Bny*uy + Bnz*uz; // B^i u_i = W B^i v_i = alpha b^0

   const double bx = (Bnlowx + ux * udotB) * invw; // b_i = P_ij B^j / w = (gamma_ij + u_i u_j) B^j / w
   const double by = (Bnlowy + uy * udotB) * invw;
   const double bz = (Bnlowz + uz * udotB) * invw;

   const double B2 = Bnlowx*Bnx + Bnlowy*Bny + Bnlowz*Bnz;
   const double b2 = (B2 + SQR(udotB)) * invw2;    // b^2 = (B^2 + (alpha b^0)^2)/W^2 from equ. 11 of arXiv:gr-qc/0701109v2 
	       
   /* Get pressure, its derivatives, and the specific enthalpy */
  double rho, press, dpdeps, dpdrho, h, dhdeps;
  rho = dens*invw;
  if ( do_polytrope == 0 ) {
       press = EOS_Pressure(handle,rho,eps);
       dpdeps = EOS_DPressByDEps(handle,rho,eps);
       dpdrho = EOS_DPressByDRho(handle,rho,eps);
       h = 1 + eps + press/rho;
       dhdeps = 1 + dpdeps/rho;
  } else {
       press = EOS_Pressure(handle,rho,1.0);
       dpdrho = EOS_DPressByDRho(handle,rho,1.0);
       dpdeps = -1e237; /* only used for do_polytrope == 0, makes gcc happy */
       h = 1 + EOS_SpecificIntEnergy(handle,rho,press) + press/rho;
       dhdeps = -1e237; /* only used for do_polytrope == 0, makes gcc happy */
  }

   /************************** 
    ** Function Evaluations ** 
    **************************/

   if ( NULL != f ) {
      /* Construct s_i definitions */
      const double fx = (dens*h + w*b2)*ux - udotB*bx - sx;
      const double fy = (dens*h + w*b2)*uy - udotB*by - sy;
      const double fz = (dens*h + w*b2)*uz - udotB*bz - sz;
      gsl_vector_set (f, 0, fx);
      gsl_vector_set (f, 1, fy);
      gsl_vector_set (f, 2, fz);

      /* Construct tau definition if not polytrope */
      if ( do_polytrope == 0 ) {
           const double fs  = dens*(w*h-1) + b2*SQR(w) - ( press + b2/2.) - SQR(udotB) - tau;
           gsl_vector_set (f, 3, fs);
      }
   }

   /************************** 
    ** Jacobian Evaluations ** 
    **************************/

   if ( NULL != J ) {
      /* Shorthands to make eqns readable */
      double dwdux = velx;
      double dwduy = vely;
      double dwduz = velz;
      double dhdux = ( press/rho - dpdrho ) * dwdux * invw;
      double dhduy = ( press/rho - dpdrho ) * dwduy * invw;
      double dhduz = ( press/rho - dpdrho ) * dwduz * invw;
      double db2dux = 2.*udotB*Bnx*invw2 - 2.*b2*dwdux*invw;
      double db2duy = 2.*udotB*Bny*invw2 - 2.*b2*dwduy*invw;
      double db2duz = 2.*udotB*Bnz*invw2 - 2.*b2*dwduz*invw;
      double dbxdux = ( ux*Bnx - bx*dwdux + udotB)*invw;
      double dbyduy = ( uy*Bny - by*dwduy + udotB)*invw;
      double dbzduz = ( uz*Bnz - bz*dwduz + udotB)*invw;
      double dbxduy = ( ux*Bny - bx*dwduy)*invw;
      double dbxduz = ( ux*Bnz - bx*dwduz)*invw; 
      double dbydux = ( uy*Bnx - by*dwdux)*invw;
      double dbyduz = ( uy*Bnz - by*dwduz)*invw;
      double dbzdux = ( uz*Bnx - bz*dwdux)*invw;
      double dbzduy = ( uz*Bny - bz*dwduy)*invw;


      double dfxdux = (dens*h + b2*w) + dens*ux*dhdux + w*ux*db2dux + b2*ux*dwdux
              - bx*Bnx - udotB*dbxdux;
      double dfyduy = (dens*h + b2*w) + dens*uy*dhduy + w*uy*db2duy + b2*uy*dwduy
              - by*Bny - udotB*dbyduy;
      double dfzduz = (dens*h + b2*w) + dens*uz*dhduz + w*uz*db2duz + b2*uz*dwduz
              - bz*Bnz - udotB*dbzduz;

      double dfxduy = dens*ux*dhduy + w*ux*db2duy + b2*ux*dwduy - bx*Bny - udotB*dbxduy;
      double dfxduz = dens*ux*dhduz + w*ux*db2duz + b2*ux*dwduz - bx*Bnz - udotB*dbxduz;
      double dfydux = dens*uy*dhdux + w*uy*db2dux + b2*uy*dwdux - by*Bnx - udotB*dbydux;
      double dfyduz = dens*uy*dhduz + w*uy*db2duz + b2*uy*dwduz - by*Bnz - udotB*dbyduz;
      double dfzdux = dens*uz*dhdux + w*uz*db2dux + b2*uz*dwdux - bz*Bnx - udotB*dbzdux;
      double dfzduy = dens*uz*dhduy + w*uz*db2duy + b2*uz*dwduy - bz*Bny - udotB*dbzduy;

      gsl_matrix_set (J, 0, 0, dfxdux);
      gsl_matrix_set (J, 0, 1, dfxduy);
      gsl_matrix_set (J, 0, 2, dfxduz);
      gsl_matrix_set (J, 1, 0, dfydux);
      gsl_matrix_set (J, 1, 1, dfyduy);
      gsl_matrix_set (J, 1, 2, dfyduz);
      gsl_matrix_set (J, 2, 0, dfzdux);
      gsl_matrix_set (J, 2, 1, dfzduy);
      gsl_matrix_set (J, 2, 2, dfzduz);

      /* epsilon derivatives if not polytrope */
      if ( do_polytrope == 0 ) {
           double dfxdeps = dens*ux*dhdeps;
           double dfydeps = dens*uy*dhdeps;
           double dfzdeps = dens*uz*dhdeps;
           double dfsdeps = dens*w*dhdeps - dpdeps;
           gsl_matrix_set (J, 0, 3, dfxdeps);
           gsl_matrix_set (J, 1, 3, dfydeps);
           gsl_matrix_set (J, 2, 3, dfzdeps);
           gsl_matrix_set (J, 3, 3, dfsdeps);

           double dfsdux = dens*dwdux*(h + press/rho + ((1.-SQR(w))*invw2)*dpdrho) 
           	+ 2*w*b2*dwdux + (SQR(w)-0.5)*db2dux - 2*udotB*Bnx;
           double dfsduy = dens*dwduy*(h + press/rho + ((1.-SQR(w))*invw2)*dpdrho) 
           	+ 2*w*b2*dwduy + (SQR(w)-0.5)*db2duy - 2*udotB*Bny;
           double dfsduz = dens*dwduz*(h + press/rho + ((1.-SQR(w))*invw2)*dpdrho) 
           	+ 2*w*b2*dwduz + (SQR(w)-0.5)*db2duz - 2*udotB*Bnz;
           gsl_matrix_set (J, 3, 0, dfsdux);
           gsl_matrix_set (J, 3, 1, dfsduy);
           gsl_matrix_set (J, 3, 2, dfsduz);
      }
   }

   return GSL_SUCCESS;

}

static int MHDsys_getUE_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
   return MHDsys_getUE_fdf(x, p, NULL, J);
}
