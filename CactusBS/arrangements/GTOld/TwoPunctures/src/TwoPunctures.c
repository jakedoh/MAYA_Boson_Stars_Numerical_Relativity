// TwoPunctures:  File  "TwoPunctures.c"
// $Id: TwoPunctures.c,v 1.6 2004/10/25 13:41:10 herrmann Exp $

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

static inline double pow2 (const double x)
{
  return x*x;
}

static inline double pow4 (const double x)
{
  return x*x*x*x;
}

static int firstcall=1;

// -------------------------------------------------------------------
void
TwoPunctures (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  typedef enum GRID_SETUP_METHOD gsm_type;
  gsm_type gsm;

  int nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int i, j, k, ntotal = n1 * n2 * n3 * nvar;
  static double *F = NULL;
  static derivs u, v;
  int precomputed=0;

  CCTK_REAL x_grid_moved;

  if (firstcall) {
    *par_b_pcmp=-1;
    *par_m_plus_pcmp=-1;
    *par_m_minus_pcmp=-1;
    *par_P_plus_pcmp1=-1;
    *par_P_plus_pcmp2=-1;
    *par_P_plus_pcmp3=-1;
    *par_P_minus_pcmp1=-1;
    *par_P_minus_pcmp2=-1;
    *par_P_minus_pcmp3=-1;
    *par_S_plus_pcmp1=-1;
    *par_S_plus_pcmp2=-1;
    *par_S_plus_pcmp3=-1;
    *par_S_minus_pcmp1=-1;
    *par_S_minus_pcmp2=-1;
    *par_S_minus_pcmp3=-1;

    firstcall=0;
  }

  if (CCTK_EQUALS(grid_setup_method, "Taylor expansion"))
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (CCTK_EQUALS(grid_setup_method, "evaluation"))
  {
    gsm = GSM_evaluation;
  }
  else
  {
    CCTK_WARN (0, "internal error");
  }

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,"precomputed values: m+=%f m-=%f",
                                 *par_m_plus_pcmp,*par_m_minus_pcmp);
    CCTK_VInfo(CCTK_THORNSTRING,"              need: m+=%f m-=%f",
                                 par_m_plus,par_m_minus);
  }

  precomputed=0;
  if (*par_b_pcmp==par_b &&
      *par_m_plus_pcmp==par_m_plus &&
      *par_m_minus_pcmp==par_m_minus &&
      *par_P_plus_pcmp1==par_P_plus[0] &&
      *par_P_plus_pcmp2==par_P_plus[1] &&
      *par_P_plus_pcmp3==par_P_plus[2] &&
      *par_P_minus_pcmp1==par_P_minus[0] &&
      *par_P_minus_pcmp2==par_P_minus[1] &&
      *par_P_minus_pcmp3==par_P_minus[2] &&
      *par_S_plus_pcmp1==par_S_plus[0] &&
      *par_S_plus_pcmp2==par_S_plus[1] &&
      *par_S_plus_pcmp3==par_S_plus[2] &&
      *par_S_minus_pcmp1==par_S_minus[0] &&
      *par_S_minus_pcmp2==par_S_minus[1] &&
      *par_S_minus_pcmp3==par_S_minus[2] ) {

      precomputed=1;
    }

  if ( F == NULL || precomputed==0 ) {
    /* Solve only when called for the first time and exactly this set
     * of parameters has already been computed */
    /* free old values if they exists */
    if(F != NULL) // assume u,v,w are allocated at the same time!
    {
	free_dvector (F, 0, ntotal - 1);
        free_derivs (&u, ntotal);
	free_derivs (&v, ntotal);
    }

    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);

    CCTK_INFO ("Solving puncture equation");
    Newton (nvar, n1, n2, n3, v, Newton_tol, Newton_maxit, cctk_time);
    F_of_v (nvar, n1, n2, n3, v, F, u, cctk_time); 

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    *admMass = (par_m_plus + par_m_minus
               - 4*par_b*PunctEvalAtArbitPosition(v.d0, 1, 0, 0, n1, n2, n3));
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass is %g", *admMass);

    /* Evaluate the ADM mass at radius given by parameter admrad */
    CCTK_REAL masstmp,radtmp;
    masstmp =PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, admrad, 0.0, 0.0);
    masstmp+=PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, 0.0, admrad, 0.0);
    masstmp+=PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, 0.0, 0.0, admrad);
    radtmp=1/sqrt(3.0)*admrad;
    masstmp+=PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, radtmp, radtmp, radtmp);
    masstmp/=4.0;

    *admMassRad = par_m_plus + par_m_minus+2.0*admrad*masstmp;
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass from r= %.19g is %.19g", admrad,*admMassRad);

    /* ADM mass for m_plus/m_minus (evaluated on puncture) */
    double u_plus,u_minus; 
    switch (gsm)
    {
      case GSM_Taylor_expansion:
        u_plus = PunctTaylorExpandAtArbitPosition
           (0, nvar, n1, n2, n3, v, par_b, 0.0, 0.0);
        u_minus = PunctTaylorExpandAtArbitPosition
           (0, nvar, n1, n2, n3, v, -par_b, 0.0, 0.0);
        break;
      case GSM_evaluation:
        u_plus = PunctIntPolAtArbitPosition
            (0, nvar, n1, n2, n3, v, par_b, 0.0, 0.0);
        u_minus = PunctIntPolAtArbitPosition
            (0, nvar, n1, n2, n3, v, -par_b, 0.0, 0.0);
        break;
      default:
        fprintf(stderr,"we don't have a known gsm: gsm=%d\n",gsm);
        fprintf(stderr,"  and taylor=%d and eval=%d\n",GSM_Taylor_expansion,GSM_evaluation);
        assert (0);
    }

    *admMass_plus = (1 + u_plus )*par_m_plus +par_m_plus*par_m_minus/(4*par_b);
    *admMass_minus =(1 + u_minus)*par_m_minus+par_m_plus*par_m_minus/(4*par_b);
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass for puncture: m+=%g",*admMass_plus);
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass for puncture: m-=%g",*admMass_minus);

    /* store the paramter values we just computed the ID */
    *par_b_pcmp=par_b;
    *par_m_plus_pcmp=par_m_plus;
    *par_m_minus_pcmp=par_m_minus;
    *par_P_plus_pcmp1=par_P_plus[0];
    *par_P_plus_pcmp2=par_P_plus[1];
    *par_P_plus_pcmp3=par_P_plus[2];
    *par_P_minus_pcmp1=par_P_minus[0];
    *par_P_minus_pcmp2=par_P_minus[1];
    *par_P_minus_pcmp3=par_P_minus[2];
    *par_S_plus_pcmp1=par_S_plus[0];
    *par_S_plus_pcmp2=par_S_plus[1];
    *par_S_plus_pcmp3=par_S_plus[2];
    *par_S_minus_pcmp1=par_S_minus[0];
    *par_S_minus_pcmp2=par_S_minus[1];
    *par_S_minus_pcmp3=par_S_minus[2];
  }

  CCTK_INFO ("Interpolating result");
  fprintf(stderr,"moving origin in x direction by %f\n",move_origin_x);

  for (k = 0; k < cctk_lsh[2]; ++k)
  {
    for (j = 0; j < cctk_lsh[1]; ++j)
    {
      for (i = 0; i < cctk_lsh[0]; ++i)
      {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

	/* make sure we are not exactly on a puncture */
	if (x[ind]==par_b && y[ind]==0 && z[ind]==0) {
	  x_grid_moved = x[ind]+ move_origin_x +1e-15;
	}
	else {
	  x_grid_moved = x[ind]+ move_origin_x;
	}
/*        fprintf(stderr,"x[ind] %f move_origin_x %f x_grid_move %f\n",x[ind],move_origin_x,x_grid_moved);*/

        const double r_plus_tmp
          = sqrt(pow2(x_grid_moved - par_b) + pow2(y[ind]) + pow2(z[ind]));
        const double r_minus_tmp
          = sqrt(pow2(x_grid_moved + par_b) + pow2(y[ind]) + pow2(z[ind]));
        
	/* modify r_plus and r_minus to get smooth behaviour on puncture */
	const double r_plus
	  = pow(pow(r_plus_tmp,4)+pow(epsilon,4),0.25);
	const double r_minus
	  = pow(pow(r_minus_tmp,4)+pow(epsilon,4),0.25);


        double U;
        switch (gsm)
        {
        case GSM_Taylor_expansion:
          U = PunctTaylorExpandAtArbitPosition
            (0, nvar, n1, n2, n3, v, x_grid_moved, y[ind], z[ind]);
          break;
        case GSM_evaluation:
          U = PunctIntPolAtArbitPosition
            (0, nvar, n1, n2, n3, v, x_grid_moved, y[ind], z[ind]);
          break;
        default:
          assert (0);
        }

        const double psi1 = 1
          + 0.5 * par_m_plus / r_plus
          + 0.5 * par_m_minus / r_minus + U;

        double static_psi = 1;
        
        double Aij[3][3];
        BY_Aijofxyz (x_grid_moved, y[ind], z[ind], cctk_time, Aij);
        
        if (keep_u_around)
            puncture_u[ind] = U;

        gxx[ind] = pow4 (psi1 / static_psi);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow4 (psi1 / static_psi);
        gyz[ind] = 0;
        gzz[ind] = pow4 (psi1 / static_psi);

        kxx[ind] = Aij[0][0] / pow2(psi1);
        kxy[ind] = Aij[0][1] / pow2(psi1);
        kxz[ind] = Aij[0][2] / pow2(psi1);
        kyy[ind] = Aij[1][1] / pow2(psi1);
        kyz[ind] = Aij[1][2] / pow2(psi1);
        kzz[ind] = Aij[2][2] / pow2(psi1);

        if (CCTK_Equals(initial_lapse,"utb")) {
          alp[ind] = 1./( pow2(psi1) );
        }
        else if (CCTK_Equals(initial_lapse,"2/(1+psi_{BL}^4)")) {
          alp[ind] = 2./( 1. + pow4(psi1) );
        }

      }
    }
  }
  CCTK_INFO("done with initial data setup for now");
  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    free_derivs (&u, ntotal);
    free_derivs (&v, ntotal);
  }
}
