#include <cctk.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "GRHydro_Analysis_initialize.h"

/* -----------------------------------------------

   M A G N E T I C    D I A G N O S T I C 
   C O M P U T A T I O N 

----------------------------------------------  */

void GRHydroAnalysisComputeMagneticDiagnostic_local(CCTK_ARGUMENTS);
void GRHydroAnalysisComputeMagneticDiagnostic_global(CCTK_ARGUMENTS);


void GRHydroAnalysisComputeMagneticDiagnostic_local(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if ( Verbosity > 0 ) {
       char strbuffer[256];
       sprintf(strbuffer,"Compute local MagneticDiagnostic at level: dx is %g (CU)",CCTK_DELTA_SPACE(0));
       CCTK_INFO(strbuffer);
    } 


//    CCTK_REAL * restrict M0   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 0) ] ); 
//    CCTK_REAL * restrict M    = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 1) ] ); 
//    CCTK_REAL * restrict J    = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 2) ] ); 
//    CCTK_REAL * restrict T    = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 3) ] ); 
//    CCTK_REAL * restrict En   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 4) ] ); 

    CCTK_REAL * restrict Emag      = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 0) ] );
    CCTK_REAL * restrict EmagTOR   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 1) ] );
    CCTK_REAL * restrict EmagPOL   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 2) ] );
    CCTK_REAL * restrict BtorNORM1 = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 3) ] );
    CCTK_REAL * restrict BpolNORM1 = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 4) ] );
    CCTK_REAL * restrict BtorM0re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 5) ] );
    CCTK_REAL * restrict BtorM1re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 6) ] );
    CCTK_REAL * restrict BtorM1im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 7) ] );
    CCTK_REAL * restrict BtorM2re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 8) ] );
    CCTK_REAL * restrict BtorM2im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 9) ] );
    CCTK_REAL * restrict BtorM3re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,10) ] );
    CCTK_REAL * restrict BtorM3im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,11) ] );
    CCTK_REAL * restrict BtorM4re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,12) ] );
    CCTK_REAL * restrict BtorM4im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,13) ] );
    CCTK_REAL * restrict BpolM0re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,14) ] );
    CCTK_REAL * restrict BpolM1re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,15) ] );
    CCTK_REAL * restrict BpolM1im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,16) ] );
    CCTK_REAL * restrict BpolM2re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,17) ] );
    CCTK_REAL * restrict BpolM2im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,18) ] );
    CCTK_REAL * restrict BpolM3re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,19) ] );
    CCTK_REAL * restrict BpolM3im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,20) ] );
    CCTK_REAL * restrict BpolM4re  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,21) ] );
    CCTK_REAL * restrict BpolM4im  = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,22) ] );


    CCTK_REAL * restrict velx = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,0)]);
    CCTK_REAL * restrict vely = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,1)]);
    CCTK_REAL * restrict velz = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,2)]);

    CCTK_REAL * restrict Bx = &(Bvec[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,0)]);
    CCTK_REAL * restrict By = &(Bvec[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,1)]);
    CCTK_REAL * restrict Bz = &(Bvec[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,2)]);

    int i,j,k;
#pragma omp for 
    for(i=0; i<cctk_lsh[0]; i++)
      for(j=0; j<cctk_lsh[1]; j++)
	for(k=0; k<cctk_lsh[2]; k++)
	  {
            CCTK_INT i3D  = CCTK_GFINDEX3D(cctkGH, i, j, k);
            // CCTK_INT i3Dx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,0);
            // CCTK_INT i3Dy = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,1);
            // CCTK_INT i3Dz = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,2);
	    
	    CCTK_REAL LOCgxx,LOCgxy,LOCgxz,LOCgyy,LOCgyz,LOCgzz;
            CCTK_REAL LOCalp,LOCbetax,LOCbetay,LOCbetaz;
            CCTK_REAL LOCrho,LOCeps,LOCpress,LOCvelx,LOCvely,LOCvelz;
            CCTK_REAL LOCx,LOCy,LOCz,LOCRcyl;

            CCTK_REAL LOCsqrtg,LOCW,LOCV2;
            CCTK_REAL LOCu0,LOCux,LOCuy,LOCuz,LOCuD0,LOCuDx,LOCuDy,LOCuDz;
            CCTK_REAL LOCbetaDx,LOCbetaDy,LOCbetaDz;
            CCTK_REAL LOCrhoStar,LOCe0;
            CCTK_REAL TMPRcyl,TMPsin1,TMPcos1,TMPsin2,TMPcos2,TMPsin3,TMPcos3,TMPsin4,TMPcos4;

	    CCTK_REAL LOCBx,LOCBy,LOCBz,LOCBDx,LOCBDy,LOCBDz;        // B field contra-variant and covariant components
	    CCTK_REAL LOCB2,LOCbst,LOCb2,LOCEMt00,LOCEmag;           // Quantities for computing the Magnetic Energy
	    CCTK_REAL LOCBphi,LOCBDphi;                              // Toroidal component
	    CCTK_REAL LOCBtorNORM2,LOCBpolNORM2,LOCBtorNORM1,LOCBpolNORM1; // Toroidal and Poloidal quantities
	    CCTK_REAL LOCBtor,LOCBpolx,LOCBpoly,LOCBpolz,LOCBpol2;   // B field components in the 'comoving splitting'
            CCTK_REAL LOCEMt00pol,LOCEMt00tor,LOCEmagPOL,LOCEmagTOR; // Magnetic energy in the 'comoving splitting'


            LOCx = x[i3D];
            LOCy = y[i3D];
            LOCz = z[i3D];

            LOCalp = alp[i3D];
            LOCbetax = betax[i3D];
            LOCbetay = betay[i3D];
            LOCbetaz = betaz[i3D];
            LOCgxx = gxx[i3D];
            LOCgxy = gxy[i3D];
            LOCgxz = gxz[i3D];
            LOCgyy = gyy[i3D];
            LOCgyz = gyz[i3D];
            LOCgzz = gzz[i3D];


	    LOCBx = Bx[i3D];
	    LOCBy = By[i3D];
	    LOCBz = Bz[i3D];

            LOCrho  = rho[i3D];
            LOCeps  = eps[i3D];
            LOCpress= press[i3D];
            LOCvelx = velx[i3D];
            LOCvely = vely[i3D];
            LOCvelz = velz[i3D];
 
            // ---------------------------------------------
            // 
            // ---------------------------------------------

            LOCsqrtg = sqrt(LOCgxx * LOCgyy * LOCgzz
               + 2.0 * LOCgxy * LOCgxz * LOCgyz
               - LOCgxy * LOCgxy * LOCgzz
               - LOCgxz * LOCgxz * LOCgyy 
               - LOCgyz * LOCgyz * LOCgxx);

           LOCV2 = (  LOCgxx * LOCvelx * LOCvelx 
		     + LOCgyy * LOCvely * LOCvely
		     + LOCgzz * LOCvelz * LOCvelz
		     + 2.0 * LOCgxy * LOCvelx * LOCvely
		     + 2.0 * LOCgxz * LOCvelx * LOCvelz
		     + 2.0 * LOCgyz * LOCvely * LOCvelz  
		    );
            LOCW = 1.0 /sqrt(1 - LOCV2);

            //  -------------------------------------------
            //  The contra-variant components of the 4-velocity
            //  :>> u^\mu
            //  ------------------------------------------
            LOCu0 = LOCW / LOCalp;
            LOCux = LOCW * (LOCvelx - LOCbetax / LOCalp);
            LOCuy = LOCW * (LOCvely - LOCbetay / LOCalp);
            LOCuz = LOCW * (LOCvelz - LOCbetaz / LOCalp);
            //  -------------------------------------------
            //  The covariant components of the 4-velocity
            //  :>> u_\mu
            //  ------------------------------------------
            LOCuD0 = LOCW * ( -LOCalp 
               + LOCgxx * LOCvelx * LOCbetax + LOCgxy * LOCvelx * LOCbetay + LOCgxz * LOCvelx * LOCbetaz 
               + LOCgxy * LOCvely * LOCbetax + LOCgyy * LOCvely * LOCbetay + LOCgyz * LOCvely * LOCbetaz 
               + LOCgxz * LOCvelz * LOCbetax + LOCgyz * LOCvelz * LOCbetay + LOCgzz * LOCvelz * LOCbetaz 
               );
            LOCuDx = LOCW * (LOCgxx * LOCvelx + LOCgxy * LOCvely + LOCgxz * LOCvelz);
            LOCuDy = LOCW * (LOCgxy * LOCvelx + LOCgyy * LOCvely + LOCgyz * LOCvelz);
            LOCuDz = LOCW * (LOCgxz * LOCvelx + LOCgyz * LOCvely + LOCgzz * LOCvelz);
            //  -------------------------------------------
            //  The covariant components of the spatial
            //  shift :>> \beta_i
            //  ------------------------------------------
            LOCbetaDx = LOCgxx * LOCbetax + LOCgxy * LOCbetay + LOCgxz * LOCbetaz;
            LOCbetaDy = LOCgxy * LOCbetax + LOCgyy * LOCbetay + LOCgyz * LOCbetaz;
            LOCbetaDz = LOCgxz * LOCbetax + LOCgyz * LOCbetay + LOCgzz * LOCbetaz;
            //  -------------------------------------
            //  Now we defined the conserved density
            //  and some of the component of mater 
            //  energy momentum tensor
            //  -------------------------------------
            
            LOCrhoStar = LOCsqrtg * LOCW * LOCrho;
            LOCe0      = LOCrho * (1 + LOCeps) + LOCpress ;

            //  Plese not that e0 is \rho*h (h=hentalpy)
            // -----------------------------------------------------
            //  Real computation of the variables to be integrated
            // -----------------------------------------------------
            // Now we are ready to define quantity useful
            // to define the modes and circular components
            // -------------------------------------------
            // // sin# = sin(#*\phi)
            // // cos# = cos(#*\phi)
            // // -------------------------------------------

            TMPRcyl = sqrt(LOCx * LOCx + LOCy * LOCy);
            TMPsin1 = ( (TMPRcyl > GRHydro_Analysis_FuzzyTiny) ? LOCy / TMPRcyl : 0.0 ) ;
            TMPcos1 = ( (TMPRcyl > GRHydro_Analysis_FuzzyTiny) ? LOCx / TMPRcyl : 0.0 ) ;
            TMPsin2 = 2 * TMPcos1 * TMPsin1;                    
            TMPcos2 = TMPcos1 * TMPcos1 - TMPsin1 * TMPsin1;  
            TMPsin3 = TMPcos1 * TMPsin2 + TMPsin1 * TMPcos2;
            TMPcos3 = TMPcos1 * TMPcos2 - TMPsin1 * TMPsin2;
            TMPsin4 = 2.0 * TMPcos2 * TMPsin2;
            TMPcos4 = TMPcos2 * TMPcos2 - TMPsin2 * TMPsin2;

	    //  -------------------------------------------
            //  The covariant components of the B field
            //  :>> B_i
            //  ------------------------------------------
	    LOCBDx = LOCgxx * LOCBx + LOCgxy * LOCBy + LOCgxz * LOCBz;
            LOCBDy = LOCgxy * LOCBx + LOCgyy * LOCBy + LOCgyz * LOCBz;
            LOCBDz = LOCgxz * LOCBx + LOCgyz * LOCBy + LOCgzz * LOCBz;

	    //  -------------------------------------------
	    //  Computations for the Magnetic Energy ...
	    //  -------------------------------------------
	    LOCB2  = LOCBx * LOCBDx + LOCBy * LOCBDy + LOCBz * LOCBDz;
	    LOCbst = LOCW * (LOCBDx * LOCvelx + LOCBDy * LOCvely + LOCBDz * LOCvelz);
	    LOCb2  = (LOCB2 + LOCbst * LOCbst) / (LOCW * LOCW);
	    LOCEMt00 = LOCB2 - 0.5 * LOCb2;

	    //  -------------------------------------------
	    //  Toroidal and Poloidal components ...
	    //  -------------------------------------------
	    LOCBphi  = ( (TMPRcyl > GRHydro_Analysis_FuzzyTiny) ? 
			 (LOCx * LOCBy - LOCy * LOCBx) / (TMPRcyl*TMPRcyl) 
			 : 0.0 ) ; 
	    LOCBDphi = LOCx * LOCBDy - LOCy * LOCBDx;
	    LOCBtorNORM2 = LOCBphi * LOCBDphi;
	    LOCBpolNORM2 = LOCB2 - LOCBtorNORM2;

	    //  -------------------------------------------
	    // A proposal for a toroidal-poloidal decomposition
	    // in the non-axisymmetric case ("comoving splitting")
	    //  -------------------------------------------
      	    LOCBtor  = ( (LOCV2 > GRHydro_Analysis_FuzzyTiny) ?
                         (LOCBDx * LOCvelx + LOCBDy * LOCvely + LOCBDz * LOCvelz) / sqrt(LOCV2) 
	    		 : 0.0 ) ;
	    LOCBpolx = LOCBx - ( (LOCV2 > GRHydro_Analysis_FuzzyTiny) ?
                                 ( LOCBtor * LOCvelx / sqrt(LOCV2) )
                                 : 0.0 ) ;
	    LOCBpoly = LOCBy - ( (LOCV2 > GRHydro_Analysis_FuzzyTiny) ?
			         ( LOCBtor * LOCvely / sqrt(LOCV2) )
			         : 0.0 ) ;
	    LOCBpolz = LOCBz - ( (LOCV2 > GRHydro_Analysis_FuzzyTiny) ?
			         ( LOCBtor * LOCvelz / sqrt(LOCV2) )
			         : 0.0 ) ;
	    LOCBpol2 = (  LOCgxx * LOCBpolx * LOCBpolx
			+ LOCgyy * LOCBpoly * LOCBpoly
			+ LOCgzz * LOCBpolz * LOCBpolz
			+ 2.0 * LOCgxy * LOCBpolx * LOCBpoly
		        + 2.0 * LOCgxz * LOCBpolx * LOCBpolz
		        + 2.0 * LOCgyz * LOCBpoly * LOCBpolz
		       );
	    LOCEMt00pol = 0.5 * (1 + LOCV2) * LOCBpol2;
	    LOCEMt00tor = 0.5 * LOCBtor * LOCBtor;


            // -----------------------------------------------------------
            // The following quantities are local scalar and, indeed, 
            // they integated  values is effected by the impossed
            // symmetry only through the corrected UnitVolume
            // -----------------------------------------------------------

	    Emag [i3D] = LOCEMt00 * LOCsqrtg;
	    EmagTOR [i3D] = LOCEMt00tor * LOCsqrtg;
	    EmagPOL [i3D] = LOCEMt00pol * LOCsqrtg;
	    BtorNORM1 [i3D] = sqrt( fabs(LOCBtorNORM2) );
	    BpolNORM1 [i3D] = sqrt( fabs(LOCBpolNORM2) );




            // -----------------------------------------------------------
            // The following quantities are quite involved and the
            // Symmetries correction is better performed at this level
            // -----------------------------------------------------------
	    if ( GlobalGRHydroAnalysisDetect.symRot180 == 1 ) {
	        BtorM0re [i3D] = fabs(LOCBtor);
	        BtorM1re [i3D] = 0.0;
	        BtorM1im [i3D] = 0.0;
	        BtorM2re [i3D] = fabs(LOCBtor) * TMPcos2 ;
	        BtorM2im [i3D] = fabs(LOCBtor) * TMPsin2 ;
	        BtorM3re [i3D] = 0.0;
	        BtorM3im [i3D] = 0.0;
	        BtorM4re [i3D] = fabs(LOCBtor) * TMPcos4 ;
	        BtorM4im [i3D] = fabs(LOCBtor) * TMPsin4 ;
	        BpolM0re [i3D] = sqrt(LOCBpol2);
		BpolM1re [i3D] = 0.0;
	        BpolM1im [i3D] = 0.0;
	        BpolM2re [i3D] = sqrt(LOCBpol2) * TMPcos2 ;
	        BpolM2im [i3D] = sqrt(LOCBpol2) * TMPsin2 ;
	        BpolM3re [i3D] = 0.0;
	        BpolM3im [i3D] = 0.0;
	        BpolM4re [i3D] = sqrt(LOCBpol2) * TMPcos4 ;
	        BpolM4im [i3D] = sqrt(LOCBpol2) * TMPsin4 ;
            } 
	    else if ( GlobalGRHydroAnalysisDetect.symRot90  == 1 ) {
	        BtorM0re [i3D] = fabs(LOCBtor);
	        BtorM1re [i3D] = 0.0;
	        BtorM1im [i3D] = 0.0; 
	        BtorM2re [i3D] = 0.0;
	        BtorM2im [i3D] = 0.0;
	        BtorM3re [i3D] = 0.0;
	        BtorM3im [i3D] = 0.0;
	        BtorM4re [i3D] = fabs(LOCBtor) * TMPcos4 ;
	        BtorM4im [i3D] = fabs(LOCBtor) * TMPsin4 ;
	        BpolM0re [i3D] = sqrt(LOCBpol2);
		BpolM1re [i3D] = 0.0;
	        BpolM1im [i3D] = 0.0;
	        BpolM2re [i3D] = 0.0;
	        BpolM2im [i3D] = 0.0;
	        BpolM3re [i3D] = 0.0;
	        BpolM3im [i3D] = 0.0;
	        BpolM4re [i3D] = sqrt(LOCBpol2) * TMPcos4 ;
	        BpolM4im [i3D] = sqrt(LOCBpol2) * TMPsin4 ;


	    } 
	    else if ( GlobalGRHydroAnalysisDetect.symXref   == 1 ) {
	        BtorM0re [i3D] = fabs(LOCBtor);
	        BtorM1re [i3D] = 0.0 ;
	        BtorM1im [i3D] = fabs(LOCBtor) * TMPsin1 ;
	        BtorM2re [i3D] = fabs(LOCBtor) * TMPcos2 ;
	        BtorM2im [i3D] = 0.0 ;
                BtorM3re [i3D] = 0.0 ;
                BtorM3im [i3D] = fabs(LOCBtor) * TMPsin3 ;
                BtorM4re [i3D] = fabs(LOCBtor) * TMPcos4 ;
                BtorM4im [i3D] = 0.0 ;
		BpolM0re [i3D] = sqrt(LOCBpol2);
		BpolM1re [i3D] = 0.0 ;
                BpolM1im [i3D] = sqrt(LOCBpol2) * TMPsin1 ;
                BpolM2re [i3D] = sqrt(LOCBpol2) * TMPcos2 ;
                BpolM2im [i3D] = 0.0 ;
                BpolM3re [i3D] = 0.0 ;
                BpolM3im [i3D] = sqrt(LOCBpol2) * TMPsin3 ;
                BpolM4re [i3D] = sqrt(LOCBpol2) * TMPcos4 ;
                BpolM4im [i3D] = 0.0 ;


	    } 

	    else if ( GlobalGRHydroAnalysisDetect.symYref   == 1 ) {
	        BtorM0re [i3D] = fabs(LOCBtor);
	        BtorM1re [i3D] = fabs(LOCBtor) * TMPcos1 ;
                BtorM1im [i3D] = 0.0 ;
                BtorM2re [i3D] = fabs(LOCBtor) * TMPcos2 ;
                BtorM2im [i3D] = 0.0 ;
                BtorM3re [i3D] = fabs(LOCBtor) * TMPcos3 ;
                BtorM3im [i3D] = 0.0 ; 
                BtorM4re [i3D] = fabs(LOCBtor) * TMPcos4 ;
                BtorM4im [i3D] = 0.0 ; 
		BpolM0re [i3D] = sqrt(LOCBpol2);
		BpolM1re [i3D] = sqrt(LOCBpol2) * TMPcos1 ;
                BpolM1im [i3D] = 0.0 ;
                BpolM2re [i3D] = sqrt(LOCBpol2) * TMPcos2 ;
                BpolM2im [i3D] = 0.0 ;
                BpolM3re [i3D] = sqrt(LOCBpol2) * TMPcos3 ;
                BpolM3im [i3D] = 0.0 ; 
                BpolM4re [i3D] = sqrt(LOCBpol2) * TMPcos4 ;
                BpolM4im [i3D] = 0.0 ; 

	    } 
	    
	    else {
	       BtorM0re [i3D] = fabs(LOCBtor);
	        BtorM1re [i3D] = fabs(LOCBtor) * TMPcos1 ;
                BtorM1im [i3D] = fabs(LOCBtor) * TMPsin1 ;
                BtorM2re [i3D] = fabs(LOCBtor) * TMPcos2 ;
                BtorM2im [i3D] = fabs(LOCBtor) * TMPsin2 ;
                BtorM3re [i3D] = fabs(LOCBtor) * TMPcos3 ;
                BtorM3im [i3D] = fabs(LOCBtor) * TMPsin3 ;
                BtorM4re [i3D] = fabs(LOCBtor) * TMPcos4 ;
                BtorM4im [i3D] = fabs(LOCBtor) * TMPsin4 ;
		BpolM0re [i3D] = sqrt(LOCBpol2);
		BpolM1re [i3D] = sqrt(LOCBpol2) * TMPcos1 ;
                BpolM1im [i3D] = sqrt(LOCBpol2) * TMPsin1 ;
                BpolM2re [i3D] = sqrt(LOCBpol2) * TMPcos2 ;
                BpolM2im [i3D] = sqrt(LOCBpol2) * TMPsin2 ;
                BpolM3re [i3D] = sqrt(LOCBpol2) * TMPcos3 ;
                BpolM3im [i3D] = sqrt(LOCBpol2) * TMPsin3 ;
                BpolM4re [i3D] = sqrt(LOCBpol2) * TMPcos4 ;
                BpolM4im [i3D] = sqrt(LOCBpol2) * TMPsin4 ;


	    }
	  }
}

void GRHydroAnalysisComputeMagneticDiagnostic_global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
      
    CCTK_INT reduction_handle_sum,reduction_handle_max;

    if ( Verbosity > 0 ) {
       char strbuffer[256];
       sprintf(strbuffer,"Reduce MagneticDiagnostic_global: Corser Volume factor is %g (CU)",GlobalGRHydroAnalysisDetect.UnitVolume);
       CCTK_INFO(strbuffer);
    }  

    reduction_handle_sum = CCTK_ReductionHandle("sum");
    if (reduction_handle_sum < 0)
      CCTK_WARN(0, "Unable to get reduction handle.");

 if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_Emag     , 1,CCTK_VarIndex("TempPool::TempPoolGF1[0]")))  CCTK_WARN(0, "Error while reducing 0");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_EmagTOR  , 1,CCTK_VarIndex("TempPool::TempPoolGF1[1]")))  CCTK_WARN(0, "Error while reducing 1");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_EmagPOL  , 1,CCTK_VarIndex("TempPool::TempPoolGF1[2]")))  CCTK_WARN(0, "Error while reducing 2");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorNORM1, 1,CCTK_VarIndex("TempPool::TempPoolGF1[3]")))  CCTK_WARN(0, "Error while reducing 3");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolNORM1, 1,CCTK_VarIndex("TempPool::TempPoolGF1[4]")))  CCTK_WARN(0, "Error while reducing 4");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM0_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[5]")))  CCTK_WARN(0, "Error while reducing 5");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM1_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[6]")))  CCTK_WARN(0, "Error while reducing 6");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM1_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[7]")))  CCTK_WARN(0, "Error while reducing 7");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM2_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[8]")))  CCTK_WARN(0, "Error while reducing 8");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM2_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[9]")))  CCTK_WARN(0, "Error while reducing 9");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM3_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[10]"))) CCTK_WARN(0, "Error while reducing 10");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM3_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[11]"))) CCTK_WARN(0, "Error while reducing 11");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM4_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[12]"))) CCTK_WARN(0, "Error while reducing 12");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BtorM4_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[13]"))) CCTK_WARN(0, "Error while reducing 13");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM0_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[14]"))) CCTK_WARN(0, "Error while reducing 14");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM1_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[15]"))) CCTK_WARN(0, "Error while reducing 15");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM1_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[16]"))) CCTK_WARN(0, "Error while reducing 16");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM2_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[17]"))) CCTK_WARN(0, "Error while reducing 17");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM2_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[18]"))) CCTK_WARN(0, "Error while reducing 18");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM3_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[19]"))) CCTK_WARN(0, "Error while reducing 19");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM3_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[20]"))) CCTK_WARN(0, "Error while reducing 20");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM4_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[21]"))) CCTK_WARN(0, "Error while reducing 21");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle_sum, 1,CCTK_VARIABLE_REAL,magnetic_BpolM4_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[22]"))) CCTK_WARN(0, "Error while reducing 22"); 

    *magnetic_Emag       *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_EmagTOR    *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_EmagPOL    *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorNORM1  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolNORM1  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM0_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM1_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM1_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM2_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM2_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM3_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM3_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM4_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BtorM4_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM0_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM1_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM1_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM2_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM2_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM3_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM3_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM4_re  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *magnetic_BpolM4_im  *= GlobalGRHydroAnalysisDetect.UnitVolume;

    if ( GlobalGRHydroAnalysisDetect.symXref   == 1 ) {
    }
    if ( GlobalGRHydroAnalysisDetect.symYref   == 1 ) {
    }
    if ( GlobalGRHydroAnalysisDetect.symZref   == 1 ) {
    } 
    if ( GlobalGRHydroAnalysisDetect.symRot90  == 1 ) {
    }
    if ( GlobalGRHydroAnalysisDetect.symRot180 == 1 ) {
    }

  // Computation of the maximum of the toroidal and poloidal component - NOT IMPLEMENTED YET
    //
     reduction_handle_max = CCTK_ReductionHandle("maximum");
     if (reduction_handle_max < 0)
       CCTK_WARN(0, "Unable to get reduction handle 'maximum'.");
    
     if (CCTK_Reduce(cctkGH, -1, reduction_handle_max, 1,CCTK_VARIABLE_REAL,magnetic_Btor ,1,CCTK_VarIndex("TempPool::TempPoolGF1[5]")))  CCTK_WARN(0, "Error while reducing 5-max");
     if (CCTK_Reduce(cctkGH, -1, reduction_handle_max, 1,CCTK_VARIABLE_REAL,magnetic_Bpol ,1,CCTK_VarIndex("TempPool::TempPoolGF1[14]")))  CCTK_WARN(0, "Error while reducing 14-max");
    



}

