#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>

#define SQR(x) ((x)*(x))

void QuasiGRMHD_ResetVelocity(CCTK_ARGUMENTS);

void QuasiGRMHD_ResetVelocity(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  //DECLARE_CCTK_FUNCTIONS;

  CCTK_INT is_polytype;

  if(CCTK_EQUALS(whisky_eos_type, "Polytype"))
    is_polytype = 1;
  else if(CCTK_EQUALS(whisky_eos_type, "General"))
    is_polytype = 0;
  else
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Unsupported EOS type '%s'.", whisky_eos_type);
    return;
  }

#pragma omp parallel for
  //for(int k = cctk_nghostzones[2] ; k < cctk_lsh[2] -cctk_nghostzones[2] ; k++)
  for(int k = 0 ; k < cctk_lsh[2] ; k++)
  {
    // shorthands to address the B and E fields
    const CCTK_REAL *B[3] = {Bff1, Bff2, Bff3};
    const CCTK_REAL *E[3] = {Eff1, Eff2, Eff3};

    //for(int j = cctk_nghostzones[1] ; j < cctk_lsh[1] - cctk_nghostzones[1] ; j++)
    for(int j = 0 ; j < cctk_lsh[1] ; j++)
    {
      //for(int i = cctk_nghostzones[0] ; i < cctk_lsh[0] - cctk_nghostzones[0]; i++)
      for(int i = 0 ; i < cctk_lsh[0] ; i++)
      {
        CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL detg, B2, E2; // det(gamma_ij), B^i B_i, E^i E_i
        CCTK_REAL Blow[3];      // B_i
        CCTK_REAL Elow[3];      // E_i
        CCTK_REAL vlow[3];      // v_i
	CCTK_REAL factor;	// (1 + em_lorentz/wlorentz)^-1
        CCTK_REAL wlorentz;     // Lorentz factor
	CCTK_REAL em_lorentz;   // EM factor = sqrt( B^2 / B^2-E^2 )
        CCTK_REAL v_par[3];     // (u^i-u_perp^i)/W
        CCTK_REAL v_perp[3];    // field perp. velocity, seee astro-ph/0601410 Eq.18

        // TODO: add check to only apply if EM fields dominate
        // if(B2_over_rho < B2_over_rho_threshold)
        //   continue;

        // compute some shorthands to make the equations managable
	detg = 2*gxy[ind]*gxz[ind]*gyz[ind] + gzz[ind]*(gxx[ind]*gyy[ind] 
             - SQR(gxy[ind])) - gyy[ind]*SQR(gxz[ind]) - gxx[ind]*SQR(gyz[ind]);
	
        Blow[0] = gxx[ind]*B[0][ind] + gxy[ind]*B[1][ind] + gxz[ind]*B[2][ind];
        Blow[1] = gxy[ind]*B[0][ind] + gyy[ind]*B[1][ind] + gyz[ind]*B[2][ind];
        Blow[2] = gxz[ind]*B[0][ind] + gyz[ind]*B[1][ind] + gzz[ind]*B[2][ind];

        Elow[0] = gxx[ind]*E[0][ind] + gxy[ind]*E[1][ind] + gxz[ind]*E[2][ind];
        Elow[1] = gxy[ind]*E[0][ind] + gyy[ind]*E[1][ind] + gyz[ind]*E[2][ind];
        Elow[2] = gxz[ind]*E[0][ind] + gyz[ind]*E[1][ind] + gzz[ind]*E[2][ind];

        vlow[0] = gxx[ind]*velx[ind] + gxy[ind]*vely[ind] + gxz[ind]*velz[ind];
        vlow[1] = gxy[ind]*velx[ind] + gyy[ind]*vely[ind] + gyz[ind]*velz[ind];
        vlow[2] = gxz[ind]*velx[ind] + gyz[ind]*vely[ind] + gzz[ind]*velz[ind];

        B2 = Blow[0]*B[0][ind] + Blow[1]*B[1][ind] + Blow[2]*B[2][ind];

        E2 = Elow[0]*E[0][ind] + Elow[1]*E[1][ind] + Elow[2]*E[2][ind];

        wlorentz = 1. / sqrt(1. - (
                    vlow[0]*velx[ind] + vlow[1]*vely[ind] + vlow[2]*velz[ind]));

	// calculate parallel (to E,B) component of the velocity
	v_par[0] = v_par[1] = v_par[2] = 0;

	if ( E2 != 0 )
	{
	     v_par[0]  += E[0][ind]*E[0][ind]/E2 * vlow[0]
			+ E[0][ind]*E[1][ind]/E2 * vlow[1]
			+ E[0][ind]*E[2][ind]/E2 * vlow[2];
             v_par[1]  += E[1][ind]*E[0][ind]/E2 * vlow[0]
                        + E[1][ind]*E[1][ind]/E2 * vlow[1]
                        + E[1][ind]*E[2][ind]/E2 * vlow[2];
             v_par[2]  += E[2][ind]*E[0][ind]/E2 * vlow[0]
                        + E[2][ind]*E[1][ind]/E2 * vlow[1]
                        + E[2][ind]*E[2][ind]/E2 * vlow[2];
	}
	if ( B2 != 0 )
        {
             v_par[0]  += B[0][ind]*B[0][ind]/B2 * vlow[0]
                        + B[0][ind]*B[1][ind]/B2 * vlow[1]
                        + B[0][ind]*B[2][ind]/B2 * vlow[2];
             v_par[1]  += B[1][ind]*B[0][ind]/B2 * vlow[0]
                        + B[1][ind]*B[1][ind]/B2 * vlow[1]
                        + B[1][ind]*B[2][ind]/B2 * vlow[2];
             v_par[2]  += B[2][ind]*B[0][ind]/B2 * vlow[0]
                        + B[2][ind]*B[1][ind]/B2 * vlow[1]
                        + B[2][ind]*B[2][ind]/B2 * vlow[2];
        }

/*        v_par[0] = (E[0][ind]*E[0][ind]/E2 + B[0][ind]*B[0][ind]/B2) * vlow[0]
                 + (E[0][ind]*E[1][ind]/E2 + B[0][ind]*B[1][ind]/B2) * vlow[1]
                 + (E[0][ind]*E[2][ind]/E2 + B[0][ind]*B[2][ind]/B2) * vlow[2];
        v_par[1] = (E[1][ind]*E[0][ind]/E2 + B[1][ind]*B[0][ind]/B2) * vlow[0]
                 + (E[1][ind]*E[1][ind]/E2 + B[1][ind]*B[1][ind]/B2) * vlow[1]
                 + (E[1][ind]*E[2][ind]/E2 + B[1][ind]*B[2][ind]/B2) * vlow[2];
        v_par[2] = (E[2][ind]*E[0][ind]/E2 + B[2][ind]*B[0][ind]/B2) * vlow[0]
                 + (E[2][ind]*E[1][ind]/E2 + B[2][ind]*B[1][ind]/B2) * vlow[1]
                 + (E[2][ind]*E[2][ind]/E2 + B[2][ind]*B[2][ind]/B2) * vlow[2]; */

	// calculate perpindicular component
	v_perp[0] = v_perp[1] = v_perp[2] = 0;

	if ( B2 != 0 )
	{
            v_perp[0] = (Elow[1]*Blow[2]-Elow[2]*Blow[1])/(sqrt(detg)*B2); // alpha in numer. and in sqrt(-g) = alpha sqrt(gamma) cancel
            v_perp[1] = (Elow[2]*Blow[0]-Elow[0]*Blow[2])/(sqrt(detg)*B2);
            v_perp[2] = (Elow[0]*Blow[1]-Elow[1]*Blow[0])/(sqrt(detg)*B2);
	}
	// calculate factor in velocity equations
	if ( B2 - E2 <= 0 )
	{ 
	    factor = 0;
	}
	else
	{
	    em_lorentz = sqrt( B2 / (B2 - E2) );
	    factor = 1.0 / ( 1.0 + em_lorentz / wlorentz );
	}

        // TODO: think about some fade in as described in McKinney's paper
        //velx[ind] = factor*(v_par[0] - v_perp[0]) + v_perp[0];
        //vely[ind] = factor*(v_par[1] - v_perp[1]) + v_perp[1];
        //velz[ind] = factor*(v_par[2] - v_perp[2]) + v_perp[2];

        //printf("%d,%d,%d\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\n",i,j,k,B[0][ind],B[1][ind],B[2][ind], E[0][ind], E[1][ind], E[2][ind]);
	//fflush(stdout);

        /*if(activate_mhd)
        {
          if(is_polytype)
            Prim2ConMHDPoly(*whisky_eos_handle, 
                         gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind],
                         detg, &dens[ind], &sx[ind], &sy[ind], &sz[ind], &tau[ind], 
                         rho[ind], velx[ind], vely[ind], velz[ind], &eps[ind], 
                         &press[ind], &w_lorentz[ind]);
          else
            Prim2ConMHDGen(*whisky_eos_handle, 
                         gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind],
                         detg, &dens[ind], &sx[ind], &sy[ind], &sz[ind], &tau[ind], 
                         rho[ind], velx[ind], vely[ind], velz[ind], &eps[ind], 
                         &press[ind], &w_lorentz[ind]);
        }
        else
        {*/

	/*if ( i == 4 && j==30 && k==27 )
	printf("%d,%d,%d\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",
			 i,j,k, gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind],
                         detg, dens[ind], sx[ind], sy[ind], sz[ind], tau[ind], 
                         rho[ind], velx[ind], vely[ind], velz[ind], eps[ind], 
                         press[ind], w_lorentz[ind]);*/

          if(is_polytype)
            Prim2ConPoly(*whisky_eos_handle, 
                         gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind],
                         detg, &dens[ind], &sx[ind], &sy[ind], &sz[ind], &tau[ind], 
                         rho[ind], velx[ind], vely[ind], velz[ind], &eps[ind], 
                         &press[ind], &w_lorentz[ind]);
          else
            Prim2ConGen(*whisky_eos_handle, 
                         gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind],
                         detg, &dens[ind], &sx[ind], &sy[ind], &sz[ind], &tau[ind], 
                         rho[ind], velx[ind], vely[ind], velz[ind], eps[ind], 
                         &press[ind], &w_lorentz[ind]); 
       //}

       /* if ( i == 4 && j==30 && k==27 )
        printf("%d,%d,%d\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",
                         i,j,k, gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind],
                         detg, dens[ind], sx[ind], sy[ind], sz[ind], tau[ind],      
                         rho[ind], velx[ind], vely[ind], velz[ind], eps[ind], 
                         press[ind], w_lorentz[ind]);*/
      }
    }
  }
  printf("finished4000\n");
}
