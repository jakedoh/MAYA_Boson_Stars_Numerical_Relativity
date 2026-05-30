#include <assert.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

#include "hll.h"

// this mostly based on the paper by Del Zanna arXiv:0704.3206v1 describing the
// ECHO code
// some equations taken from CQG 23 (2006) S505-527 Neilsen et.al.
// some equations and the overall algorithm taken from Mignone in Jets from
// Young Stars III, ISBN 978-3-540-76967-5
// flux equations taken from Giacomazzo arXiv:0701109 which describes WhiskyMHD
// could also use for remainder but its sign conventions for the eigenvalues 
// differ from the ones currently used in here
//

/***********************************************************************/
/*********** internal routines *****************************************/
/***********************************************************************/
static inline CCTK_REAL hll_flux_add(CCTK_REAL flux_L, CCTK_REAL flux_R,
    CCTK_REAL lambda_L, CCTK_REAL lambda_R, CCTK_REAL U_L, CCTK_REAL U_R);
static void RiemannSolverHLL_Solve_Internal(const cGH * cctkGH, 
    CCTK_REAL * restrict const alp, CCTK_REAL * restrict const betax,
    CCTK_REAL * restrict const betay, CCTK_REAL * restrict const betaz,
    CCTK_REAL * restrict const gxx, CCTK_REAL * restrict const gxy,
    CCTK_REAL * restrict const gxz, CCTK_REAL * restrict const gyy,
    CCTK_REAL * restrict const gyz, CCTK_REAL * restrict const gzz,
    CCTK_REAL * restrict const rhoplus, CCTK_REAL * restrict const pressplus,
    CCTK_REAL * restrict const epsplus,
    CCTK_REAL * restrict const rhominus, CCTK_REAL * restrict const pressminus,
    CCTK_REAL * restrict const epsminus,
    CCTK_REAL * restrict const velxplus, CCTK_REAL * restrict const velyplus,
    CCTK_REAL * restrict const velzplus,
    CCTK_REAL * restrict const velxminus, CCTK_REAL * restrict const velyminus,
    CCTK_REAL * restrict const velzminus,
    CCTK_REAL * restrict const Bnxplus, 
    CCTK_REAL * restrict const Bnyplus, CCTK_REAL * restrict const Bnzplus,
    CCTK_REAL * restrict const Bnxminus, 
    CCTK_REAL * restrict const Bnyminus, CCTK_REAL * restrict const Bnzminus,
    CCTK_REAL * restrict const densplus, CCTK_REAL * restrict const sxplus,
    CCTK_REAL * restrict const syplus, CCTK_REAL * restrict const szplus,
    CCTK_REAL * restrict const tauplus, CCTK_REAL * restrict const psiplus,
    CCTK_REAL * restrict const densminus, CCTK_REAL * restrict const sxminus,
    CCTK_REAL * restrict const syminus, CCTK_REAL * restrict const szminus,
    CCTK_REAL * restrict const tauminus, CCTK_REAL * restrict const psiminus,
    CCTK_REAL * restrict densflux, CCTK_REAL * restrict sxflux,
    CCTK_REAL * restrict syflux, CCTK_REAL * restrict szflux,
    CCTK_REAL * restrict tauflux,
    CCTK_REAL * restrict Bnxflux, CCTK_REAL * restrict Bnyflux, 
    CCTK_REAL * restrict Bnzflux, CCTK_REAL * restrict psiflux,
    CCTK_INT * restrict const space_mask, CCTK_INT trivial, CCTK_INT type_bits,
    CCTK_INT eos_handle, 
    CCTK_INT xoff, CCTK_INT yoff, CCTK_INT zoff, 
    CCTK_INT const imin[3], CCTK_INT const imax[3]);

/***********************************************************************/
/*********** scheduled routines ****************************************/
/***********************************************************************/
void RiemannSolverHLL_Solve(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT type_bits, trivial;
  CCTK_INT fd, xoff, yoff, zoff;
  CCTK_REAL *beta[3] = {betax, betay, betaz};
  CCTK_REAL *velplus[3] = {velxplus, velyplus, velzplus};
  CCTK_REAL *velminus[3] = {velxminus, velyminus, velzminus};
  CCTK_REAL *splus[3] = {sxplus, syplus, szplus};
  CCTK_REAL *sminus[3] = {sxminus, syminus, szminus};
  CCTK_REAL *Bn[3] = {Bnx, Bny, Bnz};
  CCTK_REAL *Bnplus[3] = {Bnxplus, Bnyplus, Bnzplus};
  CCTK_REAL *Bnminus[3] = {Bnxminus, Bnyminus, Bnzminus};
  CCTK_REAL *g1[3] = {gxx, gxy, gxz};
  CCTK_REAL *g2[3] = {gxy, gyy, gyz};
  CCTK_REAL *g3[3] = {gxz, gyz, gzz};
  CCTK_REAL **g[3] = {g1, g2, g3};
  CCTK_REAL *Sflux[3] = {sxflux, syflux, szflux};
  CCTK_REAL *Bnflux[3] = {Bnxflux, Bnyflux, Bnzflux};
  CCTK_REAL *psifluxL=NULL;
  CCTK_REAL *psiplusL=NULL;
  CCTK_REAL *psiminusL=NULL;
  CCTK_INT imin[3], imax[3];

  if (clean_divergence) {
     psifluxL=divclean_psiflux;
     psiplusL=divclean_psiplus;
     psiminusL=divclean_psiminus;
  }

  // set up copies for the direction dependent quantities
  fd = *flux_direction-1; // rotating of coordinates are even permutations of xyz
  xoff = *xoffset;
  yoff = *yoffset;
  zoff = *zoffset;
  if(fd == 0 || fd == 1 || fd == 2)
  {
    const char *masks[] = {"Hydro_RiemannProblemX","Hydro_RiemannProblemY","Hydro_RiemannProblemZ"};
    CCTK_REAL **vec[] = {beta, velplus, velminus, Bn, Bnplus, Bnminus, g1, g2, g3, Sflux, Bnflux, splus, sminus};
    CCTK_REAL ***tens[] = {g};

    type_bits = SpaceMask_GetTypeBits(masks[fd]);
    trivial = SpaceMask_GetStateBits(masks[fd], "trivial");

    // transform all vectors (and tensors) by rotating (with wraparound) the
    // components the required number of times to rotate xyz into yzx or zxy
    for(size_t v = 0 ; v < sizeof(vec)/sizeof(vec[0]) ; v++)
    {
      for(int rr = 0 ; rr < fd ; rr++)
      {
        void *tmp = vec[v][0];
        vec[v][0] = vec[v][1];
        vec[v][1] = vec[v][2];
        vec[v][2] = tmp;
      }
    }
    // really the same as for vectors, just pacify the compiler with proper types
    for(size_t t = 0 ; t < sizeof(tens)/sizeof(tens[0]) ; t++)
    {
      for(int rr = 0 ; rr < fd ; rr++)
      {
        void *tmp = tens[t][0];
        tens[t][0] = tens[t][1];
        tens[t][1] = tens[t][2];
        tens[t][2] = tmp;
      }
    }

    // do work
    for(int d = 0 ; d < 3 ; d++)
    {
      imin[d] = whisky_stencil-1;
      imax[d] = cctk_lsh[d] - whisky_stencil - 1;
    }
    RiemannSolverHLL_Solve_Internal(cctkGH, alp, beta[0], beta[1], beta[2], 
        g[0][0], g[0][1], g[0][2], g[1][1], g[1][2], g[2][2],
        rhoplus, pressplus, epsplus, 
        rhominus, pressminus, epsminus, 
        velplus[0], velplus[1], velplus[2],
        velminus[0], velminus[1], velminus[2],
        Bnplus[0], Bnplus[1], Bnplus[2], Bnminus[0], Bnminus[1], Bnminus[2],
        densplus, splus[0], splus[1], splus[2], tauplus, psiplusL, 
        densminus, sminus[0], sminus[1], sminus[2], tauminus, psiminusL,
        densflux, Sflux[0], Sflux[1], Sflux[2], tauflux,
        Bnflux[0], Bnflux[1], Bnflux[2], psifluxL,
        space_mask, trivial, type_bits,
        *whisky_eos_handle, 
        xoff, yoff, zoff, imin, imax);
  }
  else
  {
    CCTK_WARN(0, "Flux direction not x,y,z");
    return; /* NOTREACHED */
  }

}

/***********************************************************************/
/*********** internal routines *****************************************/
/***********************************************************************/

// the HLL flux averaging function, Eq. 56 of Del Zanna
static inline CCTK_REAL hll_flux_add(CCTK_REAL flux_L, CCTK_REAL flux_R,
    CCTK_REAL lambda_L, CCTK_REAL lambda_R, CCTK_REAL U_L, CCTK_REAL U_R)
{
  CCTK_REAL flux;
  flux = (lambda_R*flux_L - lambda_L*flux_R + lambda_R*lambda_L*(U_R-U_L))/(lambda_R - lambda_L);
  return flux;
}

// actual Riemann solver routine. Implements Equ. 56 of Del Zanna
static void RiemannSolverHLL_Solve_Internal(const cGH * cctkGH, 
    CCTK_REAL * restrict const alp, CCTK_REAL * restrict const betax,
    CCTK_REAL * restrict const betay, CCTK_REAL * restrict const betaz,
    CCTK_REAL * restrict const gxx, CCTK_REAL * restrict const gxy,
    CCTK_REAL * restrict const gxz, CCTK_REAL * restrict const gyy,
    CCTK_REAL * restrict const gyz, CCTK_REAL * restrict const gzz,
    CCTK_REAL * restrict const rhoplus, CCTK_REAL * restrict const pressplus,
    CCTK_REAL * restrict const epsplus,
    CCTK_REAL * restrict const rhominus, CCTK_REAL * restrict const pressminus,
    CCTK_REAL * restrict const epsminus,
    CCTK_REAL * restrict const velxplus, CCTK_REAL * restrict const velyplus,
    CCTK_REAL * restrict const velzplus,
    CCTK_REAL * restrict const velxminus, CCTK_REAL * restrict const velyminus,
    CCTK_REAL * restrict const velzminus,
    CCTK_REAL * restrict const Bnxplus,
    CCTK_REAL * restrict const Bnyplus, CCTK_REAL * restrict const Bnzplus,
    CCTK_REAL * restrict const Bnxminus,
    CCTK_REAL * restrict const Bnyminus, CCTK_REAL * restrict const Bnzminus,
    CCTK_REAL * restrict const densplus, CCTK_REAL * restrict const sxplus,
    CCTK_REAL * restrict const syplus, CCTK_REAL * restrict const szplus,
    CCTK_REAL * restrict const tauplus, CCTK_REAL * restrict const psiplus,
    CCTK_REAL * restrict const densminus, CCTK_REAL * restrict const sxminus,
    CCTK_REAL * restrict const syminus, CCTK_REAL * restrict const szminus,
    CCTK_REAL * restrict const tauminus, CCTK_REAL * restrict const psiminus,
    CCTK_REAL * restrict densflux, CCTK_REAL * restrict sxflux,
    CCTK_REAL * restrict syflux, CCTK_REAL * restrict szflux,
    CCTK_REAL * restrict tauflux,
    CCTK_REAL * restrict Bnxflux, CCTK_REAL * restrict Bnyflux, 
    CCTK_REAL * restrict Bnzflux, CCTK_REAL * restrict psiflux,
    CCTK_INT * restrict const space_mask, 
    CCTK_INT trivial, CCTK_INT type_bits,
    CCTK_INT eos_handle, 
    CCTK_INT xoff, CCTK_INT yoff, CCTK_INT zoff, 
    CCTK_INT const imin[3], CCTK_INT const imax[3])
{
  /* NOTE: do not declare CCTK_ARGUMETNS here since we use fake ones that are rotated */
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k = imin[2]; k <= imax[2] ; k++)
  {
    CCTK_INT index, indexp;
    CCTK_REAL lambda_L, lambda_R;
    CCTK_REAL avg_gxx, avg_gxy, avg_gxz, avg_gyy, avg_gyz, avg_gzz;
    CCTK_REAL avg_uxx, avg_alp;
    CCTK_REAL avg_betax, avg_betay, avg_betaz;
    CCTK_REAL avg_detg, sqrt_detg;

    for(int j = imin[1] ; j <= imax[1] ; j++)
    {
      for(int i = imin[0]  ; i <= imax[0] ; i++)
      {
        // get local copies of grid variables
	index  =  CCTK_GFINDEX3D(cctkGH,i,j,k);
	indexp =  CCTK_GFINDEX3D(cctkGH,i+xoff,j+yoff,k+zoff);

        // average quantities
        avg_gxx = (gxx[index] + gxx[indexp])/2;
        avg_gxy = (gxy[index] + gxy[indexp])/2;
        avg_gxz = (gxz[index] + gxz[indexp])/2;
        avg_gyy = (gyy[index] + gyy[indexp])/2;
        avg_gyz = (gyz[index] + gyz[indexp])/2;
        avg_gzz = (gzz[index] + gzz[indexp])/2;

	avg_detg = 2*avg_gxy*avg_gxz*avg_gyz +
          avg_gzz*(avg_gxx*avg_gyy - SQR(avg_gxy)) - 
          avg_gyy*SQR(avg_gxz) - avg_gxx*SQR(avg_gyz);
        avg_uxx = (-SQR(avg_gyz) + avg_gyy*avg_gzz)/avg_detg;

        sqrt_detg = sqrt(avg_detg);

        avg_alp = (alp[index] + alp[indexp])/2;
	if ( betax == NULL ) {
		avg_betax = 0.;
		avg_betay = 0.;
		avg_betaz = 0.;
	} else {
	        avg_betax = (betax[index] + betax[indexp])/2;
	        avg_betay = (betay[index] + betay[indexp])/2;
	        avg_betaz = (betaz[index] + betaz[indexp])/2;
	}

        // If the Riemann problem is trivial, just calculate the fluxes from the 
        // left state and skip to the next cell
        if (SpaceMask_CheckStateBits(space_mask, index, type_bits, trivial))
        {
          CCTK_REAL psiplusatindex = ( psiplus != NULL ? psiplus[index] : 0);
          CCTK_REAL tmpBnxflux, tmpBnyflux, tmpBnzflux, tmppsiflux;
          hll_numerical_flux(densplus[index], sxplus[index], syplus[index],
              szplus[index], tauplus[index],
              velxplus[index], velyplus[index], velzplus[index],
              pressplus[index], avg_alp, avg_betax, avg_betay, avg_betaz, 
              avg_gxx, avg_gxy, avg_gxz, avg_gyy, avg_gyz, avg_gzz, sqrt_detg, 
              Bnxplus[index], Bnyplus[index], Bnzplus[index], psiplusatindex,
              &densflux[index], &sxflux[index], &syflux[index], &szflux[index],
              &tauflux[index], &tmpBnxflux, &tmpBnyflux, &tmpBnzflux, &tmppsiflux);
          if (psiflux != NULL)
          {
            psiflux[index] = tmppsiflux;
          } 
          if (Bnxflux != NULL)
          {
            Bnxflux[index] = tmpBnxflux;
            Bnyflux[index] = tmpBnyflux;
            Bnzflux[index] = tmpBnzflux;
          }
        } 
        else // The end of this branch is right at the bottom of the routine
        {
          // estimate for the highest speed (correct for degenerate case or no B
          // field) from Del Zanna eqs 70ff
          {
            // NOTE the naming scheme for lambda: L/R refers to Del Zanna's plus/minus
            // notation, left/right are the state (which is plus/minus otherwise in
            // Whisky)
            CCTK_REAL lambda_L_left, lambda_R_left; // characteristic speed
            CCTK_REAL lambda_L_right, lambda_R_right;
            CCTK_REAL uBnxplus, uBnxminus, uBnyplus, uBnzplus, uBnyminus, uBnzminus;

            // the speeds are written in terms of undensitized variables
            uBnxplus = Bnxplus[index]/sqrt_detg;
            uBnyplus = Bnyplus[index]/sqrt_detg;
            uBnzplus = Bnzplus[index]/sqrt_detg;
            uBnxminus = Bnxminus[indexp]/sqrt_detg;
            uBnyminus = Bnyminus[indexp]/sqrt_detg;
            uBnzminus = Bnzminus[indexp]/sqrt_detg;

            // compute wave speeds based on left and right states
            hll_lambda_prime(rhoplus[index], epsplus[index], pressplus[index], 
                velxplus[index], velyplus[index], velzplus[index],
                uBnxplus, uBnyplus, uBnzplus,
                avg_gxx, avg_gxy, avg_gxz, avg_gyy, avg_gyz, avg_gzz,
                avg_uxx, eos_handle,
                &lambda_L_left, &lambda_R_left);
            hll_lambda_prime(rhominus[indexp], epsminus[indexp], pressminus[indexp], 
                velxminus[indexp], velyminus[indexp], velzminus[indexp],
                uBnxminus, uBnyminus, uBnzminus,
                avg_gxx, avg_gxy, avg_gxz, avg_gyy, avg_gyz, avg_gzz,
                avg_uxx, eos_handle,
                &lambda_L_right, &lambda_R_right);
            
            // convert to Eulerian speed (Eq (70))
            // divide by the lapse since Whisky factors out one power of the
            // lapse in its expression for the fluxes (compare the HLLE solver)
            lambda_L_left = lambda_L_left - avg_betax / avg_alp;
            lambda_R_left = lambda_R_left - avg_betax / avg_alp;
            lambda_L_right = lambda_L_right - avg_betax / avg_alp;
            lambda_R_right = lambda_R_right - avg_betax / avg_alp;

            // find fastest speeds (remember L = minus, R = plus)
            lambda_L = min3(0,lambda_L_left, lambda_L_right);
            lambda_R = max3(0,lambda_R_left, lambda_R_right);
            // Limit if necessary
            if ( speed_limit>0 && lambda_L < -1*speed_limit ) {
#pragma omp critical
               CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                          "Whoa Nelly! Limiting left-going speed from %g to %g.",
                          lambda_L, -1.*speed_limit);
               lambda_L = -1.*speed_limit;
            }
            if ( speed_limit>0 && lambda_R > speed_limit ) {
#pragma omp critical
               CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                          "Whoa Nelly! Limiting right-going speed from %g to %g.",
                          lambda_R, 1.*speed_limit);
               lambda_R = speed_limit;
            }

          }

          // add the numerical fluxes according to Eq (56) in Del Zanna
          // (Riemann solving happens here ;-)
          {
            CCTK_REAL flux_L[9], flux_R[9];

            // use flux equation of Giacomazzo (differs by a factor of the lapse from that of Del Zanna)
            CCTK_REAL psiplusatindex = ( psiplus != NULL ? psiplus[index] : 0);
            hll_numerical_flux(densplus[index], sxplus[index], syplus[index],
                szplus[index], tauplus[index],
                velxplus[index], velyplus[index], velzplus[index],
                pressplus[index], avg_alp, avg_betax, avg_betay, avg_betaz,
                avg_gxx, avg_gxy, avg_gxz, avg_gyy, avg_gyz, avg_gzz, sqrt_detg, 
                Bnxplus[index], Bnyplus[index], Bnzplus[index], psiplusatindex,
                &flux_L[0], &flux_L[1], &flux_L[2], &flux_L[3],
                &flux_L[4], &flux_L[5], &flux_L[6], &flux_L[7], &flux_L[8]);
            CCTK_REAL psiminusatindexp = ( psiminus != NULL ? psiminus[indexp] : 0);
            hll_numerical_flux(densminus[indexp], sxminus[indexp], syminus[indexp],
                szminus[indexp], tauminus[indexp],
                velxminus[indexp], velyminus[indexp], velzminus[indexp],
                pressminus[indexp], avg_alp, avg_betax, avg_betay, avg_betaz,
                avg_gxx, avg_gxy, avg_gxz, avg_gyy, avg_gyz, avg_gzz, sqrt_detg, 
                Bnxminus[indexp], Bnyminus[indexp], Bnzminus[indexp], psiminusatindexp,
                &flux_R[0], &flux_R[1], &flux_R[2], &flux_R[3],
                &flux_R[4], &flux_R[5], &flux_R[6], &flux_R[7], &flux_R[8]);

            densflux[index] = hll_flux_add(flux_L[0], flux_R[0], lambda_L, lambda_R, densplus[index], densminus[indexp]);
            sxflux[index] = hll_flux_add(flux_L[1], flux_R[1], lambda_L, lambda_R, sxplus[index], sxminus[indexp]);
            syflux[index] = hll_flux_add(flux_L[2], flux_R[2], lambda_L, lambda_R, syplus[index], syminus[indexp]);
            szflux[index] = hll_flux_add(flux_L[3], flux_R[3], lambda_L, lambda_R, szplus[index], szminus[indexp]);
            tauflux[index] = hll_flux_add(flux_L[4], flux_R[4], lambda_L, lambda_R, tauplus[index], tauminus[indexp]);
            if (Bnxflux != NULL)
            {
              Bnxflux[index] = hll_flux_add(flux_L[5], flux_R[5], lambda_L, lambda_R, Bnxplus[index], Bnxminus[indexp]);
              Bnyflux[index] = hll_flux_add(flux_L[6], flux_R[6], lambda_L, lambda_R, Bnyplus[index], Bnyminus[indexp]);
              Bnzflux[index] = hll_flux_add(flux_L[7], flux_R[7], lambda_L, lambda_R, Bnzplus[index], Bnzminus[indexp]);
            }
            if (psiflux != NULL)
            {
              /* wave speed of divergence violation is +/- c apparently (See GRHydro) */
              const double c = 1.0;
              psiflux[index] = hll_flux_add(flux_L[8], flux_R[8], -c, c, psiplus[index], psiminus[indexp]);
            }
          }
        } // trivial_rp
      } // for i
    } // for j
  } // for k
}
