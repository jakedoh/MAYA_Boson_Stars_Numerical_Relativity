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
TwoPuncturesMatter (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  typedef enum GRID_SETUP_METHOD gsm_type;
  gsm_type gsm;

  int nvar_psi = 1, nvar_Aij = 4, ncomp_Aij = 6;
  int n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int ntotal_psi  = n1 * n2 * n3 * nvar_psi; 
  int ntotal_Aij = n1 * n2 * n3 * nvar_Aij;
  int ntotal_AijAij = n1 * n2 * n3 * 1;
  static double *F = NULL, *AijAij = NULL;
  static derivs u, v, w, A; /* order of Aij components is 11,12,13,22,23,33 */
  static derivs v_coeffs, w_coeffs, Aij_coeffs;
  int precomputed=0;
  int one_puncture;                                        //0 if there's 1 BH, 1 if there're 2
  int success;             /* return value of Newton solver */

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

  if (  (F == NULL || precomputed==0) && bypass_solve==0 ) {
    /* Solve only when called for the first time and this set
     * of parameters has not already been computed */
    // free old values if they exists
    if(F != NULL) // assume u,v,w are allocated at the same time!
    {
	free_dvector (F, 0, ntotal_psi - 1);
        free_derivs (&u, ntotal_psi);
	free_derivs (&v, ntotal_psi);
        if (hydro_field)
        {
	    free_dvector (AijAij, 0, ntotal_AijAij - 1);
            if (!zero_velocity)
            {
  	        free_derivs (&w, ntotal_Aij);
	        free_derivs (&A, ncomp_Aij*n1*n2*n3);
            }
        }
    }

    F = dvector (0, ntotal_psi - 1);
    allocate_derivs (&u, ntotal_psi);
    allocate_derivs (&v, ntotal_psi);
    if (hydro_field)
    {
        AijAij = dvector (0, ntotal_AijAij - 1);
        if (!zero_velocity)
        {
            allocate_derivs (&w, ntotal_Aij);
            allocate_derivs (&A, ncomp_Aij*n1*n2*n3);
        }
    }
    if (use_full_transform) {
      /* v_coeffs holds the Fourier/Chebyshev transform */
      v_coeffs.d0 = dvector(0, ntotal_psi - 1);
      if (hydro_field && !zero_velocity) {
        w_coeffs.d0 = dvector(0, ntotal_Aij - 1);
        Aij_coeffs.d0 = dvector(0, ncomp_Aij * n1 * n2 * n3 - 1);
      }
    }
 
    // assume user has good sense to make par_m_plus nonzero    //
    if (fabs(par_m_minus) <= 1e-12)
    {
      one_puncture = 1;
      //M_plus_i will be ADM mass of the system, M_minus_i=0   //
      CCTK_VInfo(CCTK_THORNSTRING, "One puncture detected");
    }
    else
      one_puncture = 0;

    if(hydro_field && !zero_velocity) {
      CCTK_INFO ("Solving puncture equation (Aij)");
      NonLinEquations = NonLinEquations_Aij;
      LinEquations = LinEquations_Aij;
      makeInitialGuess = makeInitialGuess_Aij;
      success = Newton (nvar_Aij, n1, n2, n3, w, Newton_tol, Newton_maxit, cctk_time);
      if(Newton_abort_if_unsuccessful && !success) {
        CCTK_WARN(0, "Newton solver failed to converge when solving for Aij");
      }
      CCTK_VInfo (CCTK_THORNSTRING, "Constructing Aij\n");
      GetAijMatter(nvar_Aij, n1, n2, n3, w, A);
      if(use_full_transform) {
        TP_FullTransform(w.d0, w_coeffs.d0, nvar_Aij, n1, n2, n3);
        TP_FullTransform(A.d0, Aij_coeffs.d0, ncomp_Aij, n1, n2, n3);
      }
    }
    CCTK_INFO ("Solving puncture equation (psi)");
    NonLinEquations = NonLinEquations_Psi;
    LinEquations = LinEquations_Psi;
    makeInitialGuess = makeInitialGuess_Psi;
    if(hydro_field) {
      NewtonParam = AijAij;
      SetAijAij(nvar_Aij, n1, n2, n3, w, cctk_time, AijAij);
    }
    success = Newton (nvar_psi, n1, n2, n3, v, Newton_tol, Newton_maxit, cctk_time);
    if(Newton_abort_if_unsuccessful && !success) {
      CCTK_WARN(0, "Newton solver failed to converge when solving for psi");
    }
    if(use_full_transform) {
      TP_FullTransform(v.d0, v_coeffs.d0, nvar_psi, n1, n2, n3);
    }
    //F_of_v (nvar_psi, n1, n2, n3, v, F, u, cctk_time);       //

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    *admMass = (par_m_plus + par_m_minus
               - 4*par_b*PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, 1, n1, n2, n3));
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass is %g", *admMass);

    /* Evaluate the ADM mass at radius given by parameter admrad */
    CCTK_REAL masstmp,radtmp;
    masstmp =PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, admrad, 0.0, 0.0);
    masstmp+=PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, 0.0, admrad, 0.0);
    masstmp+=PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, 0.0, 0.0, admrad);
    radtmp=1/sqrt(3.0)*admrad;
    masstmp+=PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, radtmp, radtmp, radtmp);
    masstmp/=4.0;

    *admMassRad = par_m_plus + par_m_minus+2.0*admrad*masstmp;
    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass from r= %.19g is %.19g", admrad,*admMassRad);

    /* ADM mass for m_plus/m_minus (evaluated on puncture) */
    double u_plus,u_minus;
    switch (gsm)
    {
      case GSM_Taylor_expansion:
        u_plus = PunctTaylorExpandAtArbitPosition
           (0, nvar_psi, n1, n2, n3, v, par_b, 0.0, 0.0);
        u_minus = PunctTaylorExpandAtArbitPosition
           (0, nvar_psi, n1, n2, n3, v, -par_b, 0.0, 0.0);
        break;
      case GSM_evaluation:
        u_plus = PunctIntPolAtArbitPosition
            (0, nvar_psi, n1, n2, n3, v, par_b, 0.0, 0.0);
        u_minus = PunctIntPolAtArbitPosition
            (0, nvar_psi, n1, n2, n3, v, -par_b, 0.0, 0.0);
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

  // Set standard BY puncture ID with U=0 
  if (bypass_solve) {
       CCTK_INFO ("NOT solving the constraints.");
       *admMass = par_m_plus + par_m_minus;
       *admMassRad = par_m_plus + par_m_minus;
       *admMass_plus = par_m_plus + par_m_plus*par_m_minus/(4*par_b);
       *admMass_minus = par_m_minus + par_m_plus*par_m_minus/(4*par_b);
       CCTK_VInfo (CCTK_THORNSTRING, "ADM mass is %g", *admMass);
       CCTK_VInfo (CCTK_THORNSTRING, "ADM mass for puncture: m+=%g",*admMass_plus);
       CCTK_VInfo (CCTK_THORNSTRING, "ADM mass for puncture: m-=%g",*admMass_minus);
  }


  CCTK_INFO ("Interpolating result");
  if (CCTK_EQUALS(metric_type, "static conformal")) {
    if (CCTK_EQUALS(conformal_storage, "factor")) {
      *conformal_state = 1;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs")) {
      *conformal_state = 2;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")) {
      *conformal_state = 3;
    }
  } else {
    *conformal_state = 0;
  }

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,"Moving origin in x direction by %f: "
                 "Punctures to appear at %f and %f\n",
                 move_origin_x,-par_b-move_origin_x,par_b-move_origin_x);
  }

  CCTK_REAL *scalarphia = CCTK_VarDataPtr(cctkGH, 0, scalarPhia);
  CCTK_REAL *scalarpia = CCTK_VarDataPtr(cctkGH, 0, scalarPia);
  CCTK_REAL *scalarphib = CCTK_VarDataPtr(cctkGH, 0, scalarPhib);
  CCTK_REAL *scalarpib = CCTK_VarDataPtr(cctkGH, 0, scalarPib);

#pragma omp parallel for
  for (int ind = 0; ind < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; ++ind)
  {
        CCTK_REAL x_grid_moved;

	/* make sure we are not exactly on a puncture */
	if ((x[ind]+move_origin_x==par_b || x[ind]+move_origin_x==-par_b) 
		&& y[ind]==0 && z[ind]==0) {
	  x_grid_moved = x[ind] + move_origin_x +1e-15;
	}
	else {
	  x_grid_moved = x[ind] + move_origin_x;
	}
/*        fprintf(stderr,"x[ind] %.19g move_origin_x %.19g x_grid_move %.19g\n",x[ind],move_origin_x,x_grid_moved); */

        const double r_plus_tmp
          = sqrt(pow2(x_grid_moved - par_b) + pow2(y[ind]) + pow2(z[ind]));
        const double r_minus_tmp
          = sqrt(pow2(x_grid_moved + par_b) + pow2(y[ind]) + pow2(z[ind]));
        
	/* modify r_plus and r_minus to get smooth behaviour on puncture */
	const double r_plus
	  = pow(pow(r_plus_tmp,4)+pow(epsilon,4),0.25);
	const double r_minus
	  = pow(pow(r_minus_tmp,4)+pow(epsilon,4),0.25);


	// avoid NaNs
	if(r_plus == 0. || r_minus == 0.)
	    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Trying to evaluate 1/r at focal point(s) '%s%s'.", r_plus==0.?"+":"",r_minus==0.?"-":"");

        double U;
        double Aij_matter[6] = {0.,0.,0.,0.,0.,0.}; /* order is 11,12,13,22,23,33 */
	if (bypass_solve) {
		U = 0;
	} else {

	        switch (gsm)
	        {
	        case GSM_Taylor_expansion:
	          U = PunctTaylorExpandAtArbitPosition
	            (0, nvar_psi, n1, n2, n3, v, x_grid_moved, y[ind], z[ind]);
                  if (hydro_field && !zero_velocity) {
                    if (use_finite_differencing_for_Aij) {
                      double delta[3] = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)};
                      FD_AijofUP (w, PunctTaylorExpandAtArbitPosition, x_grid_moved,y[ind],z[ind], delta, Aij_matter, n1,n2,n3);
                    } else {
                      for (int c = 0; c < ncomp_Aij; c++) {
                        Aij_matter[c] = PunctTaylorExpandAtArbitPosition
                                      (c, ncomp_Aij, n1, n2, n3, A, x_grid_moved,y[ind],z[ind]);
                      }
                    }
                  }
	          break;
	        case GSM_evaluation:
                  if (use_full_transform) {
                    U = PunctIntPolCoeffsAtArbitPosition
                      (0, nvar_psi, n1, n2, n3, v_coeffs, x_grid_moved, y[ind], z[ind]);
                  } else {
                    U = PunctIntPolAtArbitPosition
                      (0, nvar_psi, n1, n2, n3, v, x_grid_moved, y[ind], z[ind]);
                  }
                  if (hydro_field && !zero_velocity) {
                    if (use_finite_differencing_for_Aij) {
                      double delta[3] = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)};
                      if (use_full_transform) {
                        FD_AijofUP (w_coeffs, PunctIntPolCoeffsAtArbitPosition, x_grid_moved,y[ind],z[ind], delta, Aij_matter, n1,n2,n3);
                      } else {
                        FD_AijofUP (w, PunctIntPolAtArbitPosition, x_grid_moved,y[ind],z[ind], delta, Aij_matter, n1,n2,n3);
                      }
                    } else {
                      for (int c = 0; c < ncomp_Aij; c++) {
                        if (use_full_transform) {
                          Aij_matter[c] = PunctIntPolCoeffsAtArbitPosition
                                        (c, ncomp_Aij, n1, n2, n3, Aij_coeffs, x_grid_moved,y[ind],z[ind]);
                        } else {
                          Aij_matter[c] = PunctIntPolAtArbitPosition
                                        (c, ncomp_Aij, n1, n2, n3, A, x_grid_moved,y[ind],z[ind]);
                        }
                      }
                    }
                  }
	          break;
	        default:
	          assert (0);
	        }

	}

        /* Add a gaussian U to whatever U is already specified */
        if (add_cvPhi) {

          if ( one_puncture == 1 ) {

             for ( int ng = 0; ng<cv_ngauss; ng++) {
                 U += cv_amp[ng]*exp(-0.5*pow2( (r_plus-cv_r0[ng])/cv_sigma[ng]) );
             }

          } else {

             for ( int ng = 0; ng<cv_ngauss; ng++) {
                 U = U + cv_amp[ng]*exp(-0.5*pow2( (r_plus-cv_r0[ng])/cv_sigma[ng]) )
                     + cv_amp[ng]*exp(-0.5*pow2( (r_minus-cv_r0[ng])/cv_sigma[ng]) ) ;
             }

          }

	}

        const double psi1 = 1
          + 0.5 * par_m_plus / r_plus
          + 0.5 * par_m_minus / r_minus + U; 

        double static_psi = 1;
        
        double Aij[3][3];
        BY_Aijofxyz (x_grid_moved, y[ind], z[ind], cctk_time, Aij);
        
        Aij[0][0] += Aij_matter[0]; 
        Aij[0][1] += Aij_matter[1]; 
        Aij[0][2] += Aij_matter[2]; 
        Aij[1][0] += Aij_matter[1]; 
        Aij[1][1] += Aij_matter[3]; 
        Aij[1][2] += Aij_matter[4]; 
        Aij[2][0] += Aij_matter[2]; 
        Aij[2][1] += Aij_matter[4]; 
        Aij[2][2] += Aij_matter[5]; 

        if (*conformal_state > 0 || keep_psi_static) {

          double xp, yp, zp, rp, ir;
          double s1, s3, s5;
          double p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;

          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = x[ind] + move_origin_x - par_b;
          yp = y[ind];
          zp = z[ind];
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          ir = 1.0/rp;


          s1 = 0.5*par_m_plus*ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = x[ind] + move_origin_x + par_b;
          yp = y[ind];
          zp = z[ind];
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          ir = 1.0/rp;

//        fprintf(stderr,"x[ind] %g move_origin_x %g",ir,rp);

          s1 = 0.5*par_m_minus*ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          psi_static[ind]=p;

          if (*conformal_state >= 1) {
            static_psi = p;
            psi[ind] = static_psi;
// warning: this was needed at some time, but not anymore??
//            psi_p[ind] = static_psi;
//            psi_p_p[ind] = static_psi;
          }

          if (*conformal_state >= 2) {
            psix[ind] = px / static_psi;
            psiy[ind] = py / static_psi;
            psiz[ind] = pz / static_psi;
          }
          if (*conformal_state >= 3) {
            psixx[ind] = pxx / static_psi;
            psixy[ind] = pxy / static_psi;
            psixz[ind] = pxz / static_psi;
            psiyy[ind] = pyy / static_psi;
            psiyz[ind] = pyz / static_psi;
            psizz[ind] = pzz / static_psi;
          }

        } /* if conformal-state > 0 */

        if(keep_u_around)
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

        /* DANGER! Adding a trace to the extrinsic curvature! */
        if (add_cvK) {
  
           double cvK = 0;
           for ( int ng=0; ng<cv_ngauss; ng++) {
                 cvK += cv_amp[ng]*y[ind]*exp(-pow2(r_plus/cv_sigma[ng])/2.) + cv_amp[ng]*y[ind]*exp(-pow2(r_minus/cv_sigma[ng])/2.);
           }
 
           // Kij = Aij + 1/3 gij K
           kxx[ind] = kxx[ind] + gxx[ind]*cvK/3.;
           kxy[ind] = kxy[ind] + gxy[ind]*cvK/3.;
           kxz[ind] = kxz[ind] + gxz[ind]*cvK/3.;
           kyy[ind] = kyy[ind] + gyy[ind]*cvK/3.;
           kyz[ind] = kyz[ind] + gyz[ind]*cvK/3.;
           kzz[ind] = kzz[ind] + gzz[ind]*cvK/3.;
 
       } 
//      /* scalar fields variables */
//      we rescale the momentum

    CCTK_REAL ini_phia, ini_grad_phiax, ini_grad_phiay, ini_grad_phiaz, ini_pia, ini_phib, ini_grad_phibx, ini_grad_phiby, ini_grad_phibz, ini_pib, ini_VV;
    calcInitialScalarFields(x[ind], y[ind], z[ind], &ini_phia, &ini_grad_phiax, &ini_grad_phiay, &ini_grad_phiaz, &ini_pia, &ini_phib, &ini_grad_phibx, &ini_grad_phiby, &ini_grad_phibz, &ini_pib, &ini_VV);
    CCTK_REAL psi6;
    psi6=psi1*psi1*psi1*psi1*psi1*psi1;
    scalarphia[ind] = ini_phia;
    /*scalarpia[ind] = ini_pia;*/
    scalarpia[ind] = ini_pia/psi6;
    scalarphib[ind] = ini_phib;
    scalarpib[ind] = ini_pib/psi6; 
    /*scalarpib[ind] = ini_pib;*/

  }

  CCTK_INFO("done with initial data setup for now");
  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal_psi - 1);
    free_derivs (&u, ntotal_psi);
    free_derivs (&v, ntotal_psi);
  }
}
