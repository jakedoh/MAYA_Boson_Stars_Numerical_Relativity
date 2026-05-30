/* Main functions for converting conservatives
   to primitives */

#include <assert.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "EOS_Base.h"
#include "RecoverMHD.h"
#include <stdio.h>
#include <string.h>

/* External Functions */
extern "C" void RecoverMHD_Conservative2Primitive(CCTK_ARGUMENTS);
extern "C" void RecoverMHD_Conservative2PrimitiveBounds(CCTK_ARGUMENTS);
extern "C" void RecoverMHD_Check_Symmetry(CCTK_ARGUMENTS);
extern "C" void RecoverMHD_Con2PrimBoundsTracer(CCTK_ARGUMENTS);
extern "C" void RecoverMHD_CheckFailMask(CCTK_ARGUMENTS);
extern "C" void CCTK_FCALL CCTK_FNAME(RecoverMHD_Con2PrimBoundsTracer)(CCTK_ARGUMENTS);

extern "C" void RecoverMHD_Conservative2Primitive(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int nx,ny,nz, ierr;
  bool polytype=0;
  char puncture_specific_info[250];

  CCTK_INT type_bits, atmosphere;
  type_bits = SpaceMask_GetTypeBits("Hydro_Atmosphere");
  atmosphere = SpaceMask_GetStateBits("Hydro_Atmosphere", "in_atmosphere");
  //CCTK_INT type2_bits=-1;
  //CCTK_INT excised=-1;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  /* Are we using a polytropic EOS? */
  if ( (CCTK_Equals(whisky_eos_type,"Polytype")) )
	polytype = 1;
  
  
  CCTK_REAL local_min_tracer;
  if (use_min_tracer != 0) {
    local_min_tracer = min_tracer;
  } else {
    local_min_tracer = 0;
  }

  /* Minima from polytropic limit */
  if ( puncture_trick != 0 )
     *puncture_utau_min = *whisky_rho_min * EOS_SpecificIntEnergy( *whisky_polytrope_handle, *whisky_rho_min, 1.0);

  /* Use an appropriate version of whisky_tau_min */
  CCTK_REAL tau_min = ( *whisky_tau_min > 1e100 ) ? tau_min_initial : *whisky_tau_min;

  /* Create handle for solver type */
  typedef enum GSL_SOLVER_TYPE gsolve_type;
  gsolve_type stype;
  if ( CCTK_EQUALS(solver_type,"Newton") )
	stype=GST_Newton;
  else if ( CCTK_EQUALS(solver_type,"gNewton") )
	stype=GST_gNewton;
  else if ( CCTK_EQUALS(solver_type,"HybridSJ") )
	stype=GST_HybridSJ;
  else if ( CCTK_EQUALS(solver_type,"HybridJ") )
	stype=GST_HybridJ;
  else if ( CCTK_EQUALS(solver_type,"HybridS") )
	stype=GST_HybridS;
  else if ( CCTK_EQUALS(solver_type,"Hybrid") )
	stype=GST_Hybrid;
  else if ( CCTK_EQUALS(solver_type,"dNewton") )
	stype=GST_dNewton;
  else if ( CCTK_EQUALS(solver_type,"Broyden") )
	stype=GST_Broyden;
  else {
	CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Unknown solver type '%s'",solver_type);
	return; /* NOTREACHED */
  }

  /* Handle choice of system of eqns for recovery here to avoid */
  enum SYSTEM_EQN_TYPE{ HD_ONLY=0, MHD_5D_UiEps=1, MHD_2D_WP=2, MHD_1D_hW2=3, MHD_2D_hv2=4 } solve_system_method;
  if ( CCTK_Equals(system_type, "u_i eps") ) {
     solve_system_method = MHD_5D_UiEps;
  } else if (CCTK_Equals(system_type, "W press") ) {
     solve_system_method = MHD_2D_WP;
  } else if (CCTK_Equals(system_type, "h*W^2") ) {
     solve_system_method = MHD_1D_hW2;
  } else if (CCTK_Equals(system_type, "h v^2") ) {
     solve_system_method = MHD_2D_hv2;
  } else {
     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "unknown system type '%s'", system_type);
     return; /* NOTREACHED */
  }

  CCTK_LOOP3_ALL(recoverMHD_Pointwise,cctkGH,i,j,k) {
//  CCTK_LOOP3(recoverMHD_pointwise,
//    i,j,k, 0,0,0, nx,ny,nz,
//    cctk_ash[0],cctk_ash[1],cctk_ash[2]) 
//  #pragma omp parallel
//  for (int k=0; k < nz; k++) {
//    for (int j=0; j < ny; j++) {
//      for (int i=0; i < nx; i++) {

        /* do not compute if in atmosphere or in an excised region */
            // if (SpaceMask_CheckStateBits(space_mask, i, j, k, type_bits, atmosphere)) continue;
	/* do compute, to ensure consitent primitives (the metric might change after all) */

	/* Copy gridfunctions to local variables */
	CCTK_INT idx=CCTK_GFINDEX3D(cctkGH, i,j,k);
	CCTK_INT bxidx=CCTK_VECTGFINDEX3D(cctkGH, i,j,k,0);
	CCTK_INT byidx=CCTK_VECTGFINDEX3D(cctkGH, i,j,k,1);
	CCTK_INT bzidx=CCTK_VECTGFINDEX3D(cctkGH, i,j,k,2);

        CCTK_REAL gxxL, gxyL, gxzL, gyyL, gyzL, gzzL;
        gxxL=gxx[idx];
        gyyL=gyy[idx];
        gzzL=gzz[idx];
	gxyL=gxy[idx];
	gxzL=gxz[idx];
	gyzL=gyz[idx];
	CCTK_REAL carpetweightL = Whisky_CarpetWeights[idx];

	int do_impose_s2_limit = 0;
        CCTK_REAL rr_pct[number_of_punctures];
	if ( puncture_trick != 0 ) {
            sprintf(puncture_specific_info,"Radii relative to punctures:");
	    for (int pcn=0; pcn<number_of_punctures; pcn++) {
		 rr_pct[pcn] = SQR(x[idx] - center_x[pcn]) + 
                               SQR(y[idx] - center_y[pcn]) + 
                               SQR(z[idx] - center_z[pcn]);
                 char pctinfo[30];
                 sprintf(pctinfo," r_%d=%g(%g)",pcn,rr_pct[pcn],radii[pcn]);
                 strcat( puncture_specific_info, pctinfo );
		 if ( radii[pcn] >= 0 && (rr_pct[pcn] <= SQR(radii[pcn])) ) {
		        do_impose_s2_limit = 1;
                 }
	    }
	}

	/* conservative variables to local quantities */
	double densL, sxL, syL, szL, BnxL, BnyL, BnzL, tauL;
	densL = dens[idx];
	sxL = sx[idx];
	syL = sy[idx];
	szL = sz[idx];
	BnxL = Bnx[idx];
	BnyL = Bny[idx];
	BnzL = Bnz[idx];
	tauL = tau[idx];

	/* hydro primitives (initial guesses) */
	double rhoL, velxL, velyL, velzL, pressL, w_lorentzL, epsL,
		BvecxL, BvecyL, BveczL;
	rhoL = rho[idx];
	velxL = velx[idx];
	velyL = vely[idx];
	velzL = velz[idx];
	pressL = press[idx];
	w_lorentzL = w_lorentz[idx];
	BvecxL = Bvec[bxidx];
	BvecyL = Bvec[byidx];
	BveczL = Bvec[bzidx];
	epsL = eps[idx];
	
    	/* Apply Faber's tau limit if applicable at this point */
    	if ( do_impose_s2_limit ) {
    	    if ( tauL <= 0 && puncture_tau_min != -2 ) { // is tau nonsensical, do we care?
    	        if( carpetweightL == 1) {
    	            CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Negative tau encountered at a physical point: %g at (%g,%g,%g), %s",
    	                       tauL, x[idx], y[idx], z[idx], puncture_specific_info);
    	        }
    	    }

    	    if ( puncture_tau_min == -1 ) {
	    	double det = spatialdeterminant( gxxL, gxyL, gxzL, gyyL, gyzL, gzzL );
    	        tauL = sqrt(det) * *puncture_utau_min;
    	    } else {
    	        tauL = puncture_tau_min;
    	    }
    	}

	ierr = RecoverMHD_Con2Prim_pointwise( *whisky_eos_handle, *whisky_polytrope_handle, polytype, cctk_iteration, stype,
		(int) solve_system_method, &densL, &tauL, &sxL, &syL, &szL, &BnxL, &BnyL, &BnzL,
		&rhoL, &epsL, &velxL, &velyL, &velzL, &BvecxL, &BvecyL, &BveczL,
		&w_lorentzL, &pressL, x[idx], y[idx], z[idx],
		gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
		*whisky_reflevel, Whisky_CarpetWeights[idx], 
		do_impose_s2_limit, *whisky_rho_min, tau_min );
        if ( failed_con2prim_mask != NULL ) 
             failed_con2prim_mask[idx] += 1.*ierr;

	/* Convert tracers if active */
        if (evolve_tracer != 0) {
           for ( int itracer=0; itracer < number_of_tracers; itracer++ )
	   {
		CCTK_INT tracer_ind = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, itracer);
		tracer[tracer_ind] = cons_tracer[tracer_ind]/densL;

		if (use_min_tracer != 0) {
			if ( tracer[tracer_ind] <= local_min_tracer )
				tracer[tracer_ind] = local_min_tracer;
                }
           }
        }

	/* Copy back conservative variables to GFs in case Faber limit changed them */
	dens[idx] = densL;
	sx[idx] = sxL;
	sy[idx] = syL;
	sz[idx] = szL;
	Bnx[idx] = BnxL;
	Bny[idx] = BnyL;
	Bnz[idx] = BnzL;
	tau[idx] = tauL;

	/* Copy new primitives to GFs */
	rho[idx] = rhoL;
	velx[idx] = velxL;
	vely[idx] = velyL;
	velz[idx] = velzL;
	press[idx] = pressL;
	w_lorentz[idx] = w_lorentzL;
	Bvec[bxidx] = BvecxL;
	Bvec[byidx] = BvecyL;
	Bvec[bzidx] = BveczL;
	eps[idx] = epsL;

//      }
//    }
//  }
  }
  CCTK_ENDLOOP3_ALL(recoverMHD_Pointwise);

}

extern "C" void RecoverMHD_Conservative2PrimitiveBounds(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int nx,ny,nz, ierr;
  CCTK_REAL gxxL, gxyL, gxzL, gyyL, gyzL, gzzL;
  CCTK_REAL gxxR, gxyR, gxzR, gyyR, gyzR, gzzR;
  CCTK_REAL local_min_tracer, carpetweightL;
  bool polytype=0;
  char puncture_specific_info[250];
  

  /*CCTK_INT type_bits, atmosphere;
  type_bits = SpaceMask_GetTypeBits("Hydro_Atmosphere");
  atmosphere = SpaceMask_GetStateBits("Hydro_Atmosphere", "in_atmosphere");
  CCTK_INT type2_bits=-1;
  CCTK_INT excised=-1; */

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];
  
  if ( use_min_tracer != 0) {
    local_min_tracer = min_tracer;
  } else {
    local_min_tracer = 0.;
  }
  
  /* Are we using a polytropic EOS? */
  if ( (CCTK_Equals(whisky_eos_type,"Polytype")) )
	polytype = 1;
  
  /* Minima from polytropic limit */
  if ( puncture_trick != 0)
     *puncture_utau_min = *whisky_rho_min * EOS_SpecificIntEnergy( *whisky_polytrope_handle, *whisky_rho_min, 1.0);

  /* Use an appropriate version of whisky_tau_min */
  CCTK_REAL tau_min = ( *whisky_tau_min > 1e100 ) ? tau_min_initial : *whisky_tau_min;

  /* Create handle for solver type */
  typedef enum GSL_SOLVER_TYPE gsolve_type;
  gsolve_type stype;
  if ( CCTK_EQUALS(solver_type,"Newton") )
	stype=GST_Newton;
  else if ( CCTK_EQUALS(solver_type,"gNewton") )
	stype=GST_gNewton;
  else if ( CCTK_EQUALS(solver_type,"HybridSJ") )
	stype=GST_HybridSJ;
  else if ( CCTK_EQUALS(solver_type,"HybridJ") )
	stype=GST_HybridJ;
  else if ( CCTK_EQUALS(solver_type,"HybridS") )
	stype=GST_HybridS;
  else if ( CCTK_EQUALS(solver_type,"Hybrid") )
	stype=GST_Hybrid;
  else if ( CCTK_EQUALS(solver_type,"dNewton") )
	stype=GST_dNewton;
  else if ( CCTK_EQUALS(solver_type,"Broyden") )
	stype=GST_Broyden;
  else {
	CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,"Unknown solver type '%s'",solver_type);
	return; /* NOTREACHED */
  }

  /* Handle choice of system of eqns for recovery here to avoid */
  enum SYSTEM_EQN_TYPE{ HD_ONLY=0, MHD_5D_UiEps=1, MHD_2D_WP=2, MHD_1D_hW2=3, MHD_2D_hv2=4 } solve_system_method;
  if ( CCTK_Equals(system_type, "u_i eps") ) {
     solve_system_method = MHD_5D_UiEps;
  } else if (CCTK_Equals(system_type, "W press") ) {
     solve_system_method = MHD_2D_WP;
  } else if (CCTK_Equals(system_type, "h*W^2") ) {
     solve_system_method = MHD_1D_hW2;
  } else if (CCTK_Equals(system_type, "h v^2") ) {
     solve_system_method = MHD_2D_hv2;
  } else {
     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "unknown system type '%s'", system_type);
     return; /* NOTREACHED */
  }


  CCTK_LOOP3_ALL(recoverMHD_pointwise_cellbounds, cctkGH, i,j,k) {
//  CCTK_LOOP3(recoverMHD_pointwise_cellbounds,
//    i,j,k, 0,0,0, nx,ny,nz,
//    cctk_ash[0],cctk_ash[1],cctk_ash[2]) {
//  for (int k=0; k<nz; k++) {
//    for (int j=0; j<ny; j++) {
//      for (int i=0; i<nx; i++) {

	CCTK_INT idx=CCTK_GFINDEX3D(cctkGH, i,j,k);

        /* do not compute if in atmosphere */
        if ( SpaceMask_CheckState(space_mask, idx, "Hydro_Atmosphere", "in_atmosphere") ) 
		continue;
         
	CCTK_INT Lidx=CCTK_GFINDEX3D(cctkGH, i-*xoffset, j-*yoffset, k-*zoffset );
	CCTK_INT Ridx=CCTK_GFINDEX3D(cctkGH, i+*xoffset, j+*yoffset, k+*zoffset );

        gxxL = 0.5 * (gxx[idx] + gxx[Lidx]);
        gxyL = 0.5 * (gxy[idx] + gxy[Lidx]);
        gxzL = 0.5 * (gxz[idx] + gxz[Lidx]);
        gyyL = 0.5 * (gyy[idx] + gyy[Lidx]);
        gyzL = 0.5 * (gyz[idx] + gyz[Lidx]);
        gzzL = 0.5 * (gzz[idx] + gzz[Lidx]);

        gxxR = 0.5 * (gxx[idx] + gxx[Ridx]);
        gxyR = 0.5 * (gxy[idx] + gxy[Ridx]);
        gxzR = 0.5 * (gxz[idx] + gxz[Ridx]);
        gyyR = 0.5 * (gyy[idx] + gyy[Ridx]);
        gyzR = 0.5 * (gyz[idx] + gyz[Ridx]);
        gzzR = 0.5 * (gzz[idx] + gzz[Ridx]);

	carpetweightL = Whisky_CarpetWeights[idx];

        if (evolve_tracer != 0) {
           for ( int itracer=0; itracer < number_of_tracers; itracer++ )
	   {
		CCTK_INT tracer_ind = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, itracer);
		tracer[tracer_ind] = cons_tracer[tracer_ind] / dens[idx] ;

		if (use_min_tracer != 0) {
			if ( tracer[tracer_ind] <= local_min_tracer )
				tracer[tracer_ind] = local_min_tracer;
                }
           }
        }	
	
	/* Check whether we are to apply Faber's limit */
        /* Also, set up where we are relative to all punctures */
	int do_impose_s2_limit = 0;
        CCTK_REAL rr_pct[number_of_punctures];
	if ( puncture_trick != 0 ) {
            sprintf(puncture_specific_info,"Radii relative to punctures:");
	    for (int pcn=0; pcn<number_of_punctures; pcn++) {
		 rr_pct[pcn] = SQR(x[idx] - center_x[pcn]) + 
                               SQR(y[idx] - center_y[pcn]) + 
                               SQR(z[idx] - center_z[pcn]);
                 char pctinfo[30];
                 sprintf(pctinfo," r_%d=%g(%g)",pcn,rr_pct[pcn],radii[pcn]);
                 strcat( puncture_specific_info, pctinfo );
		 if ( radii[pcn] >= 0 && (rr_pct[pcn] <= SQR(radii[pcn])) ) {
		        do_impose_s2_limit = 1;
                 }
	    }
	}

	/*******************/
	/* Convert on left */
	/*******************/
	double densL, tauL, sxL, syL, szL, BnxL, BnyL, BnzL;
	densL = densminus[idx];
	sxL = sxminus[idx];
	syL = syminus[idx];
	szL = szminus[idx];
	BnxL = Bnxminus[idx];
	BnyL = Bnyminus[idx];
	BnzL = Bnzminus[idx];
	tauL = tauminus[idx];

	double rhoL, epsL, velxL, velyL, velzL, BvecxL, BvecyL, BveczL, w_lorentzL, pressL;
	rhoL = rhominus[idx];
	pressL = pressminus[idx];
	w_lorentzL = w_lorentzminus[idx];
	velxL = velxminus[idx];
	velyL = velyminus[idx];
	velzL = velzminus[idx];
	BvecxL = Bxminus[idx];
	BvecyL = Byminus[idx];
	BveczL = Bzminus[idx];
	epsL = epsminus[idx];

    	/* Apply Faber's limit for tau if applicable at this point */
	double avg_detL;
    	if ( do_impose_s2_limit ) {
    	    if ( tauL <= 0 && puncture_tau_min != -2 ) { // is tau nonsensical, do we care?
    	        if( carpetweightL == 1) {
    	            CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Negative tau encountered: %g at (%g,%g,%g), %s",
    	                       tauL, x[idx], y[idx], z[idx], puncture_specific_info);
    	        }
    	    }

    	    if ( puncture_tau_min == -1 ) {
	    	avg_detL = spatialdeterminant( gxxL, gxyL, gxzL, gyyL, gyzL, gzzL );
    	        tauL = sqrt(avg_detL) * *puncture_utau_min;
    	    } else {
    	        tauL = puncture_tau_min;
    	    }
    	}

	ierr = RecoverMHD_Con2Prim_pointwise( *whisky_eos_handle, *whisky_polytrope_handle, polytype, cctk_iteration, stype,
		(int) solve_system_method, &densL, &tauL, &sxL, &syL, &szL, &BnxL, &BnyL, &BnzL,
		&rhoL, &epsL, &velxL, &velyL, &velzL, &BvecxL, &BvecyL, &BveczL,
		&w_lorentzL, &pressL, x[idx] - 0.5*CCTK_DELTA_SPACE(1), 
		y[idx] - 0.5*CCTK_DELTA_SPACE(2), z[idx] - 0.5*CCTK_DELTA_SPACE(3),
		gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
		*whisky_reflevel, Whisky_CarpetWeights[idx],
		do_impose_s2_limit, *whisky_rho_min, tau_min );
        if ( failed_con2prim_mask != NULL ) 
             failed_con2prim_mask[idx] += 10.*ierr;

	/* Copy back variables to GFs */
	densminus[idx] = densL;
	sxminus[idx] = sxL;
	syminus[idx] = syL;
	szminus[idx] = szL;
	Bnxminus[idx] = BnxL;
	Bnyminus[idx] = BnyL;
	Bnzminus[idx] = BnzL;
	tauminus[idx] = tauL;

	rhominus[idx] = rhoL;
	velxminus[idx] = velxL;
	velyminus[idx] = velyL;
	velzminus[idx] = velzL;
	pressminus[idx] = pressL;
	w_lorentzminus[idx] = w_lorentzL;
	Bxminus[idx] = BvecxL;
	Byminus[idx] = BvecyL;
	Bzminus[idx] = BveczL;
	epsminus[idx] = epsL;


	/*******************/
	/* Convert on right */
	/*******************/
	double densR, tauR, sxR, syR, szR, BnxR, BnyR, BnzR;
	densR = densplus[idx];
	sxR = sxplus[idx];
	syR = syplus[idx];
	szR = szplus[idx];
	BnxR = Bnxplus[idx];
	BnyR = Bnyplus[idx];
	BnzR = Bnzplus[idx];
	tauR = tauplus[idx];

	double rhoR, epsR, velxR, velyR, velzR, BvecxR, BvecyR, BveczR, w_lorentzR, pressR;
	rhoR = rhoplus[idx];
	pressR = pressplus[idx];
	w_lorentzR = w_lorentzplus[idx];
	velxR = velxplus[idx];
	velyR = velyplus[idx];
	velzR = velzplus[idx];
	BvecxR = Bxplus[idx];
	BvecyR = Byplus[idx];
	BveczR = Bzplus[idx];
	epsR = epsplus[idx];

    	/* Apply Faber's tau limit if applicable at this point */
	double avg_detR;
    	if ( do_impose_s2_limit ) {
    	    if ( tauR <= 0 && puncture_tau_min != -2 ) { // is tau nonsensical, do we care?
    	        if( carpetweightL == 1) {
    	            CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Negative tau encountered on plus cell boundary: %g at (%g,%g,%g), %s",
    	                       tauR, x[idx], y[idx], z[idx], puncture_specific_info);
    	        }
    	    }

    	    if ( puncture_tau_min == -1 ) {
	    	avg_detR = spatialdeterminant( gxxR, gxyR, gxzR, gyyR, gyzR, gzzR );
    	        tauR = sqrt(avg_detR) * *puncture_utau_min;
    	    } else {
    	        tauR = puncture_tau_min;
    	    }
    	}

	ierr = RecoverMHD_Con2Prim_pointwise( *whisky_eos_handle, *whisky_polytrope_handle, polytype, cctk_iteration, stype,
		(int) solve_system_method, &densR, &tauR, &sxR, &syR, &szR, &BnxR, &BnyR, &BnzR,
		&rhoR, &epsR, &velxR, &velyR, &velzR, &BvecxR, &BvecyR, &BveczR,
		&w_lorentzR, &pressR, x[idx] + 0.5*CCTK_DELTA_SPACE(1), 
		y[idx] + 0.5*CCTK_DELTA_SPACE(2), z[idx] + 0.5*CCTK_DELTA_SPACE(3),
		gxxR, gxyR, gxzR, gyyR, gyzR, gzzR,
		*whisky_reflevel, Whisky_CarpetWeights[idx],
		do_impose_s2_limit, *whisky_rho_min, tau_min );

	/* Copy back variables to GFs */
	densplus[idx] = densR;
	sxplus[idx] = sxR;
	syplus[idx] = syR;
	szplus[idx] = szR;
	Bnxplus[idx] = BnxR;
	Bnyplus[idx] = BnyR;
	Bnzplus[idx] = BnzR;
	tauplus[idx] = tauR;
        if ( failed_con2prim_mask != NULL )
           failed_con2prim_mask[idx] += 100.*ierr;

	rhoplus[idx] = rhoR;
	velxplus[idx] = velxR;
	velyplus[idx] = velyR;
	velzplus[idx] = velzR;
	pressplus[idx] = pressR;
	w_lorentzplus[idx] = w_lorentzR;
	Bxplus[idx] = BvecxR;
	Byplus[idx] = BvecyR;
	Bzplus[idx] = BveczR;
	epsplus[idx] = epsR;

//      }
//    }
//  }
  }
  CCTK_ENDLOOP3_ALL(recoverMHD_pointwise_cellbounds);

} 

// wrapper routine to be called from Fortran (currenlyt only ReconstructTVD)
extern "C" void CCTK_FCALL CCTK_FNAME(RecoverMHD_Con2PrimBoundsTracer)(CCTK_ARGUMENTS)
{
  RecoverMHD_Con2PrimBoundsTracer(CCTK_PASS_CTOC);
}

extern "C" void RecoverMHD_Con2PrimBoundsTracer(CCTK_ARGUMENTS)
{ 
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int nx,ny,nz;
 
  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];
 
  CCTK_LOOP3_INT(recoverMHD_pointwise_cellboundsTracer, cctkGH, i,j,k) {
//  CCTK_LOOP3(recoverMHD_pointwise_cellboundsTracer,
//    i,j,k, whisky_stencil-1,whisky_stencil-1,whisky_stencil-1, nx - whisky_stencil,
//    ny-whisky_stencil,nz-whisky_stencil,
//    cctk_ash[0],cctk_ash[1],cctk_ash[2]) 
//  for (int k = whisky_stencil-1; k < nz - whisky_stencil; k++) {
//    for (int j = whisky_stencil-1; j < ny - whisky_stencil; j++) {
//      for (int i = whisky_stencil-1; i < nx - whisky_stencil; i++) {
        
	CCTK_INT idx=CCTK_GFINDEX3D(cctkGH, i,j,k);

	CCTK_INT Lidx=CCTK_GFINDEX3D(cctkGH, i-*xoffset, j-*yoffset, k-*zoffset );
	CCTK_INT Ridx=CCTK_GFINDEX3D(cctkGH, i+*xoffset, j+*yoffset, k+*zoffset );

        CCTK_REAL gxxL, gxyL, gxzL, gyyL, gyzL, gzzL, avg_detL;
        CCTK_REAL gxxR, gxyR, gxzR, gyyR, gyzR, gzzR, avg_detR;
        CCTK_REAL ugxxL, ugxyL, ugxzL, ugyyL, ugyzL, ugzzL;
        CCTK_REAL ugxxR, ugxyR, ugxzR, ugyyR, ugyzR, ugzzR;
        gxxL = 0.5 * (gxx[idx] + gxx[Lidx]);
        gxyL = 0.5 * (gxy[idx] + gxy[Lidx]);
        gxzL = 0.5 * (gxz[idx] + gxz[Lidx]);
        gyyL = 0.5 * (gyy[idx] + gyy[Lidx]);
        gyzL = 0.5 * (gyz[idx] + gyz[Lidx]);
        gzzL = 0.5 * (gzz[idx] + gzz[Lidx]);
        gxxR = 0.5 * (gxx[idx] + gxx[Ridx]);
        gxyR = 0.5 * (gxy[idx] + gxy[Ridx]);
        gxzR = 0.5 * (gxz[idx] + gxz[Ridx]);
        gyyR = 0.5 * (gyy[idx] + gyy[Ridx]);
        gyzR = 0.5 * (gyz[idx] + gyz[Ridx]);
        gzzR = 0.5 * (gzz[idx] + gzz[Ridx]);

	avg_detL = spatialdeterminant( gxxL, gxyL, gxzL, gyyL, gyzL, gzzL );
	avg_detR = spatialdeterminant( gxxR, gxyR, gxzR, gyyR, gyzR, gzzR );
	invertspatialdet( gxxL, gxyL, gxzL, gyyL, gyzL, gzzL, avg_detL, 
		ugxxL, ugxyL, ugxzL, ugyyL, ugyzL, ugzzL );
	invertspatialdet( gxxR, gxyR, gxzR, gyyR, gyzR, gzzR, avg_detR, 
		ugxxR, ugxyR, ugxzR, ugyyR, ugyzR, ugzzR );

	/* Check determinant */
        if ( (std::isnan(avg_detL) || avg_detL < 0) ) 
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Conservative2Primitive: Left Metric determinant in iteration %d at (%g,%g,%g) on level %d is: %g from data g=(%g,%g,%g,%g,%g,%g)", 
		cctk_iteration, x[idx],y[idx],z[idx],*whisky_reflevel,avg_detL,gxxL,gxyL,gxzL,gyyL,gyzL,gzzL);
        
        if ( (std::isnan(avg_detR) || avg_detR < 0) ) 
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Conservative2Primitive: Right Metric determinant in iteration %d at (%g,%g,%g) on level %d is: %g from data g=(%g,%g,%g,%g,%g,%g)", 
		cctk_iteration, x[idx],y[idx],z[idx],*whisky_reflevel,avg_detR,gxxR,gxyR,gxzR,gyyR,gyzR,gzzR);

	double w_lorentzR = 1./ ( sqrt(1 - ( gxxR*SQR(velxplus[idx]) + gyyR*SQR(velyplus[idx]) + gzzR*SQR(velzplus[idx]) 
		+ 2.*gxyR*velxplus[idx]*velyplus[idx] + 2.*gxzR*velxplus[idx]*velzplus[idx]
		+ 2.*gyzR*velyplus[idx]*velzplus[idx]) ) );
	double w_lorentzL = 1./ ( sqrt(1 - ( gxxL*SQR(velxminus[idx]) + gyyL*SQR(velyminus[idx]) + gzzL*SQR(velzminus[idx]) 
		+ 2.*gxyL*velxminus[idx]*velyminus[idx] + 2.*gxzL*velxminus[idx]*velzminus[idx]
		+ 2.*gyzL*velyminus[idx]*velzminus[idx]) ) );

        for ( int itracer = 0; itracer < number_of_tracers; itracer++ ) {

		/* This is just carry-over for hydro for now. Fix. */
		CCTK_INT tidx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,itracer);
		tracerplus[tidx] = cons_tracerplus[tidx] / ( sqrt(avg_detR) * w_lorentzR * rhoplus[idx] );
		tracerminus[tidx] = cons_tracerminus[tidx] / ( sqrt(avg_detL) * w_lorentzL * rhominus[idx] );

        }

/*      tracerplus[idx] = cons_tracerplus[idx] / ( sqrt(avg_detr) * rhoplus[idx] * w_lorentzR );
        tracerminus[idx] = cons_tracerminus[idx] / ( sqrt(avg_detl) * rhominus[idx] * w_lorentzL ); */
        
//      }
//    }
//  }
  }
  CCTK_ENDLOOP3_INT(recoverMHD_pointwise_cellboundsTracer);
  
}

extern "C" void RecoverMHD_Check_Symmetry(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT nx,ny,nz, gx,gy,gz, midx, pidx, cmidx, cpidx;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  gx = cctk_nghostzones[0];
  gy = cctk_nghostzones[1];
  gz = cctk_nghostzones[2];

  CCTK_INT tidx=CCTK_GFINDEX3D(cctkGH, 0,0,gz);
  if( z[tidx] > CCTK_DELTA_SPACE(3)/2. || z[tidx] < -CCTK_DELTA_SPACE(3)/2.) {
    if ( verbose > 0)
	CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Z-Symmetry check: z(gz)=%g. Skipping.",z[tidx]);
    return;
  }
  
  CCTK_LOOP3_INT(recoverMHD_checkSymm, cctkGH, i,j,k) {
//  CCTK_LOOP3(recoverMHD_checkSymm,
//    i,j,k, 0, gy, gx, gz, ny-gy, nx-gx,
//    cctk_ash[0],cctk_ash[1],cctk_ash[2]) 
//  for (int k=0; k<gz; k++) {
//    for (int j=gy; j<ny-gy; j++) {
//      for (int i=gx; i<nx-gx; i++) {

        /*if ( verbose > 1 )
           CCTK_VInfo(CCTK_THORNSTRING,"Z-symmetry comparision at %d,%d,%d. Compare %d,%d.",i,j,k,gz-k,gz+k);*/

	midx = CCTK_GFINDEX3D(cctkGH, i, j, gz-k );
	pidx = CCTK_GFINDEX3D(cctkGH, i, j, gz+k );
#define CHECK(X,Y,fact) \
        if( !check_interval(X[pidx], fact*X[midx], zsymmetry_rel_err, zsymmetry_abs_err)             \
          && Whisky_CarpetWeights[pidx]==1.) { \
          	CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Symmetry broken at index (%d,%d,%d-%d), " \
			"coordinates (%g,%g,%g) vs. (%g,%g,%g) for variable %s " \
			"on iteration %d: %g != %g.  Difference=%g", i,j,gz+k,gz-k,x[pidx],y[pidx],z[pidx], \
                        x[midx],y[midx],z[midx],Y,cctk_iteration, X[pidx],X[midx], X[pidx]-fact*X[midx]); \
        }
  
        CHECK(rho,"rho",1.0);
        CHECK(eps,"eps",1.0);
        CHECK(press,"press",1.0);
        CHECK(velx,"velx",1.0);
        CHECK(vely,"vely",1.0);
        CHECK(velz,"velz",-1.0);
        CHECK(w_lorentz,"lorentz",1.0);
        CHECK(dens,"dens",1.0);
        CHECK(tau,"tau",1.0);
        CHECK(sx,"sx",1.0);
        CHECK(sy,"sy",1.0);
        CHECK(sz,"sz",-1.0);
        CHECK(gxx,"gxx",1.0);
        CHECK(gxy,"gxy",1.0);
        CHECK(gxz,"gxz",-1.0);
        CHECK(gyy,"gyy",1.0);
        CHECK(gyz,"gyz",-1.0);
        CHECK(gzz,"gzz",1.0);
        CHECK(alp,"alp",1.0);
        CHECK(betax,"betax",1.0);
	CHECK(betay,"betay",1.0);
	CHECK(betaz,"betaz",-1.0);
        if (clean_divergence)
           CHECK(divclean_psi,"divclean_psi",1.0);

        CHECK(Bnx,"Bnx",-1.0);
        CHECK(Bny,"Bny",-1.0);
        CHECK(Bnz,"Bnz",1.0);

#define CHECK4D(X,Y,fact,vdir) \
        cpidx = CCTK_VECTGFINDEX3D(cctkGH, i, j, gz+k, vdir); \
        cmidx = CCTK_VECTGFINDEX3D(cctkGH, i, j, gz-k, vdir); \
        if( !check_interval(X[cpidx], fact*X[cmidx], zsymmetry_rel_err, zsymmetry_abs_err)             \
          && Whisky_CarpetWeights[pidx]==1.) { \
          	CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"Symmetry broken at index (%d,%d,%d/%d), " \
			"coordinates (%g,%g,%g) vs. (%g,%g,%g) for variable %s " \
			"on iteration %d: %g != %g.  Difference=%g", i,j,gz+k,gz-k,x[pidx],y[pidx],z[pidx], \
                        x[midx],y[midx],z[midx],Y,cctk_iteration, X[cpidx],X[cmidx], X[cpidx]-fact*X[cmidx]); \
        }
  
        CHECK4D(Bvec,"Bvec(x)",-1.0,0);
        CHECK4D(Bvec,"Bvec(y)",-1.0,1);
        CHECK4D(Bvec,"Bvec(z)",1.0,2);
        if ( *whisky_mhd_handle > 2 ) { 
           CHECK4D(Avec,"Avec(x)",1.0,0);
           CHECK4D(Avec,"Avec(y)",1.0,1);
           CHECK4D(Avec,"Avec(z)",-1.0,2);
        }

//      }
//    }
//  }
  }
  CCTK_ENDLOOP3_INT(recoverMHD_checkSymm);
}

extern "C" void RecoverMHD_CheckFailMask(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_INT nx, ny, nz, c2p_fail_warnlevel;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  if ( abort_on_fails > 0 ) { 
     c2p_fail_warnlevel = 0;
  } else {
     c2p_fail_warnlevel = 1;
  }

  CCTK_LOOP3_ALL(recoverMHD_checkFailMask, cctkGH, i,j,k ) {
//  CCTK_LOOP3(recoverMHD_checkFailMask,
//    i,j,k, 0, 0, 0, nz, ny, nz,
//    cctk_ash[0],cctk_ash[1],cctk_ash[2]) 
//  for ( int k = 0; k < nz; k++ )
//     for ( int j = 0; j < ny; j++ )
//        for ( int i = 0; i < nx; i++ ) {
           CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k );
           CCTK_REAL rhofac = rho[idx]/(*whisky_rho_min);

           if ( failed_con2prim_mask[idx] > 0 &&
              ( !SpaceMask_CheckState(space_mask, idx, "Hydro_Atmosphere", "in_atmosphere") &&
                Whisky_CarpetWeights[idx] > 0.5 && alp[idx] > 0.25 && rhofac > (1.+atmo_tolerance) ) ) {

              CCTK_INT x_idx = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
              CCTK_INT y_idx = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
              CCTK_INT z_idx = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);
              CCTK_REAL BxL = Bvec[x_idx];
              CCTK_REAL ByL = Bvec[y_idx];
              CCTK_REAL BzL = Bvec[z_idx];
              CCTK_REAL B2 = gxx[idx]*SQR(BxL) + gyy[idx]*SQR(ByL) + gzz[idx]*SQR(BzL)
                             + 2.*gxy[idx]*BxL*ByL + 2.*gxz[idx]*BxL*BzL + 2.*gyz[idx]*ByL*BzL;

              /* Warn. Failure codes interpreted as follows:
                    1-9    : Failure on initial recovery attempt (accumulates over substeps. 1-4. )
                    10-20  : Failure on minus side of reconstructed variables 
                    100-110: Failure on plus side of reconstructed variables */
              CCTK_VWarn(c2p_fail_warnlevel,__LINE__,__FILE__,CCTK_THORNSTRING,
                  "Con2Prim failure code %d at (%g,%g,%g) on level %d (weight %g) from conservatives "
                  " (dens,tau,sx,sy,sz,B2) = (%g,%g,%g,%g,%g,%g) with rho = %g*rho_atmo",
                  failed_con2prim_mask[idx], x[idx], y[idx], z[idx], *whisky_reflevel,
                  Whisky_CarpetWeights[idx],dens[idx], tau[idx], sx[idx], sy[idx], sz[idx], B2, rhofac );
           }

//  }
  } CCTK_ENDLOOP3_ALL(recoverMHD_checkFailMask);

}

