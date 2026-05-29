#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "EOS_Base.h"

#define SQR(x) ((x) * (x))

/* Externally Visible Functions */
void Cooler_Init(CCTK_ARGUMENTS);
void Cooler_Radiate(CCTK_ARGUMENTS);

/* Internal Functions */
static double Cool_Sharma(double pressure, double rho);

void Cooler_Init(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   for(int k = 0 ; k < cctk_lsh[2] ; k++)
   { 
     for(int j = 0 ; j < cctk_lsh[1] ; j++)
     {
       for(int i = 0 ; i < cctk_lsh[0] ; i++)
       {
         CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
         LostEps[idx] = 0.;
       }
     }
   }

}

/* Apply cooling to RHS of conservatives */
void Cooler_Radiate(CCTK_ARGUMENTS)
{

   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   int L_type;
   if ( CCTK_EQUALS(rate_function,"Sharma") ) 
	L_type=0;
   else
        CCTK_WARN(1,"Choosen Cooling rate function is not yet implemented.");


   for (int k = 0 ; k < cctk_lsh[2] ; k++)
   {
      for(int j = 0 ; j < cctk_lsh[1] ; j++)
      {
         for(int i = 0 ; i < cctk_lsh[0] ; i++)
         {
	     CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
             CCTK_REAL pressL = press[idx];
	     CCTK_REAL rhoL   = rho[idx];
	     CCTK_REAL epsL   = eps[idx];
	     CCTK_REAL wL = w_lorentz[idx];
	     CCTK_REAL vxL = velx[idx];
	     CCTK_REAL vyL = vely[idx];
	     CCTK_REAL vzL = velz[idx];

	     CCTK_REAL densL = dens[idx];

	     CCTK_REAL taurhsL = taurhs[idx];
	     CCTK_REAL sxrhsL = sxrhs[idx];
	     CCTK_REAL syrhsL = syrhs[idx];
	     CCTK_REAL szrhsL = szrhs[idx];

	     CCTK_REAL gxxL = gxx[idx];
	     CCTK_REAL gxyL = gxy[idx];
	     CCTK_REAL gxzL = gxz[idx];
	     CCTK_REAL gyyL = gyy[idx];
	     CCTK_REAL gyzL = gyz[idx];
	     CCTK_REAL gzzL = gzz[idx];
	     CCTK_REAL detg = - SQR(gxzL)*gyyL - gxxL*SQR(gyzL) - SQR(gxyL)*gzzL
			      + 2*gxyL*gxzL*gyzL + gxxL*gyyL*gzzL;

	     CCTK_REAL dtau=0;
	     CCTK_REAL dsx=0;
	     CCTK_REAL dsy=0;
	     CCTK_REAL dsz=0;

	     CCTK_REAL dPdeps = EOS_DPressByDEps(*whisky_eos_handle,rhoL,epsL);

	     CCTK_REAL rad_rate=0;
	     switch (L_type)
	     {
	        case 0:   /* Sharma */
		   rad_rate = Cool_Sharma(pressL,rhoL);
		   break;
		default:
		   CCTK_WARN(1,"Cooling rate choice error. This should have been caught earlier.");
		   return; /* Not reached */
	     }

	     dtau = densL * wL * (1. + dPdeps/rhoL) * rad_rate - sqrt(detg)*dPdeps*rad_rate;
	     dsx = densL * wL * (1. + dPdeps/rhoL) * rad_rate * (gxxL * vxL + gxyL * vyL + gxzL * vzL); 
	     dsy = densL * wL * (1. + dPdeps/rhoL) * rad_rate * (gxyL * vxL + gyyL * vyL + gyzL * vzL); 
	     dsz = densL * wL * (1. + dPdeps/rhoL) * rad_rate * (gxzL * vxL + gyzL * vyL + gzzL * vzL); 

	     if ( verbose && ( (idx%100) ==0 ) ) {
		CCTK_VInfo(CCTK_THORNSTRING,"Cooling for location (%g,%g,%g): rad_rate=%g.", x[idx],y[idx],z[idx],rad_rate);
		CCTK_VInfo(CCTK_THORNSTRING,"            dtau=%g ( %g tau when dtau/dt=%g )",dtau,dtau/tau[idx],taurhsL);
		CCTK_VInfo(CCTK_THORNSTRING,"            ds_i=(%g,%g,%g): (%g,%g,%g) s_i when ds_i/dt=(%g,%g,%g)",dsx,dsy,dsz,
			   dsx/sx[idx],dsy/sy[idx],dsz/sz[idx],sxrhsL,syrhsL,szrhsL);
	     }

	     /* Copy back to GF */
	     taurhs[idx] = taurhsL - dtau;
	     sxrhs[idx] = sxrhsL - dsx;
	     syrhs[idx] = syrhsL - dsy;
	     szrhs[idx] = szrhsL - dsz;
	     if ( keep_radiation ) {
	     	LostEps[idx] = rad_rate;
	     }

         }
      }
   }

}

/* Apply cooling to primitives to RHS of conservatives */
void Cool_Primitives(CCTK_ARGUMENTS)
{

   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   int L_type;
   if ( CCTK_EQUALS(rate_function,"Sharma") ) 
	L_type=0;
   else
        CCTK_WARN(1,"Choosen Cooling rate function is not yet implemented.");

   if ( CCTK_EQUALS(whisky_eos_type,"Polytype") )
	CCTK_WARN(0,"Cooling does not work for polytype equations of state.");

   CCTK_REAL pmin = EOS_Pressure(*whisky_polytrope_handle,*whisky_rho_min,1.);

   for (int k = 0 ; k < cctk_lsh[2] ; k++)
   {
      for(int j = 0 ; j < cctk_lsh[1] ; j++)
      {
         for(int i = 0 ; i < cctk_lsh[0] ; i++)
         {
	     CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
             CCTK_REAL pressL = press[idx];
	     CCTK_REAL rhoL   = rho[idx];
	     CCTK_REAL epsL   = eps[idx];
	     CCTK_REAL wL = w_lorentz[idx];
	     CCTK_REAL vxL = velx[idx];
	     CCTK_REAL vyL = vely[idx];
	     CCTK_REAL vzL = velz[idx];

	     CCTK_REAL densL = dens[idx];
	     CCTK_REAL tauL = tau[idx];
	     CCTK_REAL sxL = sx[idx];
	     CCTK_REAL syL = sy[idx];
	     CCTK_REAL szL = sz[idx];

	     CCTK_REAL gxxL = gxx[idx];
	     CCTK_REAL gxyL = gxy[idx];
	     CCTK_REAL gxzL = gxz[idx];
	     CCTK_REAL gyyL = gyy[idx];
	     CCTK_REAL gyzL = gyz[idx];
	     CCTK_REAL gzzL = gzz[idx];
	     CCTK_REAL detg = - SQR(gxzL)*gyyL - gxxL*SQR(gyzL) - SQR(gxyL)*gzzL
			      + 2*gxyL*gxzL*gyzL + gxxL*gyyL*gzzL;

	     CCTK_REAL dPdeps = EOS_DPressByDEps(*whisky_eos_handle,rhoL,epsL);

	     CCTK_REAL delta_eps=0;
	     switch (L_type)
	     {
	        case 0:   /* Sharma */
		   delta_eps = Cool_Sharma(pressL,rhoL);
		   break;
		default:
		   CCTK_WARN(1,"Cooling rate choice error. This should have been caught earlier.");
		   return; /* Not reached */
	     }

	     CCTK_REAL eps_old = epsL;
	     epsL = eps_old - delta_eps*CCTK_DELTA_TIME;
	     pressL = EOS_Pressure(*whisky_eos_handle,rhoL,epsL);

	     /* Make sure we're not making eps go negative */
	     if ( epsL < 0 ) {
		pressL = pmin;
		epsL = EOS_SpecificIntEnergy(*whisky_polytrope_handle,rhoL,pressL);
		delta_eps = eps_old;
	     }
 
	     if ( verbose && ( (idx%100) ==0 ) && eps_old>0 ) {
		if (eps_old>0) {
			CCTK_VInfo(CCTK_THORNSTRING,"Cooling for location (%g,%g,%g): delta_eps=%g, which is %g eps.", 
				x[idx],y[idx],z[idx],delta_eps,delta_eps/eps_old );
		} else {
			CCTK_VInfo(CCTK_THORNSTRING,"Cooling for location (%g,%g,%g): delta_eps=%g when eps=%g.", 
				x[idx],y[idx],z[idx],delta_eps,eps_old );
		}
		
	     }

	     /* Regenerate conservatives from these primitives */
	     Prim2ConGen( *whisky_eos_handle, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
			  detg, &densL, &sxL, &syL, &szL, &tauL, rhoL, vxL, 
			  vyL, vzL, epsL, &pressL, &wL );

	     /* Copy all back to GFs */
	     eps[idx] = epsL;
	     press[idx] = pressL;
	     w_lorentz[idx] = wL;
	     if ( keep_radiation ) {
	     	LostEps[idx] = delta_eps;
	     }

	     dens[idx] = densL;
	     sx[idx] = sxL;
	     sy[idx] = syL;
	     sz[idx] = szL;
	     tau[idx] = tauL;

         }
      }
   }

}

/* Implement various cooling functions */

static double Cool_Sharma(double pressure, double rho)
{

   DECLARE_CCTK_PARAMETERS;
   double lambda_cgs, Tkev, mkeV, rad_rate;

   mkeV = ion_type * 9.368e5; 		// rest energy of unit gas particle in keV
   Tkev = Tfac * mkeV * (pressure/rho);

   if ( Tkev > 0.02 ) {
	lambda_cgs = 1e-22 * ( 0.6e-3*pow(Tkev,-1.7) + 0.058*pow(Tkev,0.5)+0.063);
   } else {
	if ( Tkev < 0.0017235 )
		lambda_cgs = 1.544e-22 * pow(Tkev/0.0017235,6);
   	else
        	lambda_cgs = 6.72e-22 * pow(Tkev/0.02,0.6);
   }

   //rad_rate = 3.55e39*SQR(rho)*lambda_cgs / ( Msys_SolarMasses * SQR(ion_type) ); 
   rad_rate = 1.207e39*rho*lambda_cgs / ( Msys_SolarMasses * SQR(ion_type) ); 
   return rad_rate;

}
