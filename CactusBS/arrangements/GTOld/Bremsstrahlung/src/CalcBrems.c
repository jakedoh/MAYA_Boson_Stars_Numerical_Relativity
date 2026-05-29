#define INITVALUE (42)

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Reduction.h"
#include "HydroBase.h"

#define KRANC_C
#include "GenericFD.h"
#include "Symmetry.h"

#define SQR(x) ((x) * (x))
#define PI 3.14159265358979

//*******************************//
//****** External Routines ******//
//*******************************//
void Brems_Scalars(CCTK_ARGUMENTS);
void Brems_Init(CCTK_ARGUMENTS);
void BremsCalc(CCTK_ARGUMENTS);
void ZeroBrems(CCTK_ARGUMENTS);
void SumBrems(CCTK_ARGUMENTS);
void Brems_Output(CCTK_ARGUMENTS);

//*******************************//
//****** Internal Routines ******//
//*******************************//
static int sumLocalBrems(CCTK_ARGUMENTS, CCTK_REAL * locVolBS, 
   CCTK_REAL * locVolRBS, CCTK_REAL * locVolDBS, CCTK_REAL * locVolDRBS);
static void bremsstrahlung_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc);

static CCTK_INT file_created;
static CCTK_INT n_DTerms;

void Brems_Scalars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT   is_symbnd[6];

  /*** Initialize ***/
  *SymmetryFactor = INITVALUE;
  for (int vol=0; vol<intvolumes; vol++) {
      TotBrems[vol] = INITVALUE;
      TotRelBrems[vol] = INITVALUE;
      for ( int nobs=0; nobs<num_obs_angles; nobs++ ) {
          int ind2d=vol+intvolumes*nobs;
          TotDopBrems[ind2d] = INITVALUE;
          TotDopRelBrems[ind2d] = INITVALUE;
      }
  }
  file_created=0;
  n_DTerms=0;
  for ( int sign_el=0; sign_el<24; sign_el++ ) {
      Sym_Signs[sign_el]=100;
  }

  /*** Symmetry info: Non-doppler-boosted ***/
  int ierr = GetSymmetryBoundaries(cctkGH,6,is_symbnd);
  if ( ierr < 0 ) {
        CCTK_WARN(1,"Error in returning the symmetries.");
  }
  assert(!ierr);

  CCTK_INT symfactL=1;
  for (int l = 0; l < 6 ; l++)
  {
        if (is_symbnd[l]) symfactL = symfactL * 2;
  }

  CCTK_VInfo(CCTK_THORNSTRING,"Due to symmetries, we multiply by a factor of %d",symfactL);
  *SymmetryFactor = symfactL;

  /*** Symmetry info: Doppler-boosted ***/
  CCTK_INT Syms[3]={0,0,0};
  CCTK_INT set_signs = 0;
  for ( int d=0; d < (int)strlen(grid_symmetries); d++ ) {
      if ( grid_symmetries[d] == 'x' )
	 Syms[0] =  1;
      if ( grid_symmetries[d] == 'y' )
	 Syms[1] =  1;
      if ( grid_symmetries[d] == 'z' )
	 Syms[2] =  1;
      if ( grid_symmetries[d] == 'r' )
	 Syms[0] = -1;
  }
  n_DTerms = pow( 2, abs(Syms[0])+abs(Syms[1])+abs(Syms[2]) );
  for ( int icoord=0; icoord<3; icoord++ ) { /* First set of signs always positive. */
      Sym_Signs[icoord]=1;
      set_signs = 1;
  }

  if ( Syms[2] == 1 ) { /* Z reflections first */
     /* Copy elements present, changing relevant sign for reflection in x */
     for ( int el=0; el<(3*set_signs); el++ ) {

         if ( (el+1)%3 == 0 )
            Sym_Signs[3*set_signs+el] = -1*Sym_Signs[el];
         else
            Sym_Signs[3*set_signs+el] = Sym_Signs[el]; /* Copy */

     }
     set_signs=set_signs*2;
  }

  if ( Syms[1] == 1 ) { /* Y reflections next */
     /* Copy elements present, changing relevant sign for reflection in x */
     for ( int el=0; el<(3*set_signs); el++ ) {

         if ( (el+2)%3 == 0 )
            Sym_Signs[3*set_signs+el] = -1*Sym_Signs[el];
         else
            Sym_Signs[3*set_signs+el] = Sym_Signs[el]; /* Copy */

     }
     set_signs=set_signs*2;
  }

  if ( Syms[0] != 0 ) { 
     /* Copy elements present, changing relevant sign for reflection in x */
     for ( int el=0; el<(3*set_signs); el++ ) {

         Sym_Signs[3*set_signs+el] = Sym_Signs[el]; /* Copy */
         if ( el%3 == 0 )
            Sym_Signs[3*set_signs+el] = -1*Sym_Signs[el];
         if ( (el+2)%3 == 0 )
            Sym_Signs[3*set_signs+el] = Syms[0]*Sym_Signs[el];

     }
     set_signs=set_signs*2;
  }
  
  if ( set_signs != n_DTerms ) {
     CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "We set %d sets of signs for the Doppler shift. We expected %d sets. Syms=(%d,%d,%d)"
                 "Sym_Signs={ (%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),"
                 "(%d,%d,%d)}",set_signs, n_DTerms, Syms[0], Syms[1], Syms[2],Sym_Signs[0],Sym_Signs[1],
		 Sym_Signs[2],Sym_Signs[3],Sym_Signs[4],Sym_Signs[5],Sym_Signs[6],Sym_Signs[7],Sym_Signs[8],
		 Sym_Signs[9],Sym_Signs[10],Sym_Signs[11],Sym_Signs[12],Sym_Signs[13],Sym_Signs[14],Sym_Signs[15],
		 Sym_Signs[16],Sym_Signs[17],Sym_Signs[18],Sym_Signs[19],Sym_Signs[20],Sym_Signs[21],Sym_Signs[22],
		 Sym_Signs[23]);
  }
  if ( verbose ) {
     CCTK_VInfo ( CCTK_THORNSTRING, "n_DTerms=%d. Syms=(%d,%d,%d). "
                 "Sym_Signs={ (%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),(%d,%d,%d),"
                 "(%d,%d,%d)}", n_DTerms, Syms[0], Syms[1], Syms[2],Sym_Signs[0],Sym_Signs[1],
		 Sym_Signs[2],Sym_Signs[3],Sym_Signs[4],Sym_Signs[5],Sym_Signs[6],Sym_Signs[7],Sym_Signs[8],
		 Sym_Signs[9],Sym_Signs[10],Sym_Signs[11],Sym_Signs[12],Sym_Signs[13],Sym_Signs[14],Sym_Signs[15],
		 Sym_Signs[16],Sym_Signs[17],Sym_Signs[18],Sym_Signs[19],Sym_Signs[20],Sym_Signs[21],Sym_Signs[22],
		 Sym_Signs[23]);
  }

  /*** Volume and Emissivity Factors (code vs cgs units) ***/
  CCTK_REAL dVFact;
  if (do_cgs_units) {
      /* EmissFact includes initial 1.4e-27, sqrt(1.09e13) from T, and (6.17e17)^2 from rho conversion */
      *EmissFact=1.7596e15 * gaunt_factor * SQR(Z_gas/(SQR(Msys_SolarMasses)*gas_massunit));
      dVFact=3.24e15 * Msys_SolarMasses*SQR(Msys_SolarMasses);
  } else {
      *EmissFact=unitsfudge;
      dVFact=1;
  }
  CCTK_VInfo(CCTK_THORNSTRING,"Bremsstrahlung emissivities are multiplied by %g, dV factor is %g",*EmissFact,dVFact);

}

void Brems_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for(int k = 0 ; k < cctk_lsh[2] ; k++)
  {
    for(int j = 0 ; j < cctk_lsh[1] ; j++)
    {
      for(int i = 0 ; i < cctk_lsh[0] ; i++)
      {
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
	Brems[idx] = INITVALUE;
	RelBrems[idx] = INITVALUE;
      }
    }
  }

}

void BremsCalc(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rhoL, epsL, BremsL, RelBremsL;

  CCTK_REAL Te_13L, Tfac;

  Tfac = temperature_fac * (eos_gamma - 1.);

  for(int k=0 ; k<cctk_lsh[2] ; k++)
  {
    for(int j=0 ; j<cctk_lsh[1] ; j++)
    {
      for(int i=0 ; i<cctk_lsh[0] ; i++)
      {

        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);

	rhoL = INITVALUE;
	epsL = INITVALUE;

	BremsL = INITVALUE;
	RelBremsL = INITVALUE;

	rhoL = rho[idx];
	epsL = eps[idx];

	/* Te_13 = 1e-13 T_e(Kelvin) ~ ( kB/m_p c^2 ) T_e(Kelvin) */
	Te_13L = Tfac*epsL;

	BremsL = *EmissFact * SQR(rhoL) * sqrt( Te_13L);
	RelBremsL = BremsL * (1. + 4.796e3*Te_13L);

	Brems[idx] = BremsL;
	RelBrems[idx] = RelBremsL;
	
      }
    }
  } 

}

void ZeroBrems(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vol, nobs;

  for (vol=0;vol<intvolumes;vol++)
  {
     TotBrems[vol] = 0;
     TotRelBrems[vol] = 0;
     for ( nobs=0; nobs<num_obs_angles; nobs++ )
     {
       int ind2d=vol+intvolumes*nobs;
       TotDopBrems[ind2d] = 0;
       TotDopRelBrems[ind2d] = 0;
     }
  }
}


void SumBrems(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr, vol, nobs;
  int do_calc_bs, my_calc_bs_every;

  do_calc_bs = 0;
  for (vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_bs_every = calc_bs_every_vol[vol] >= 0 ? calc_bs_every_vol[vol] : calc_bs_every;
    assert(my_calc_bs_every >= 0);
    if ( (my_calc_bs_every > 0) && (cctk_iteration % my_calc_bs_every == 0) )
    {
        do_calc_bs = 1;
        break;
    }
  }
  if (do_calc_bs) 
  {

     double locVolBS[intvolumes];
     double locVolRBS[intvolumes];
     double locVolDBS[intvolumes*num_obs_angles];
     double locVolDRBS[intvolumes*num_obs_angles];

     for ( vol=0; vol<intvolumes; vol++)
     {
       locVolBS[vol] = -1e44;;
       locVolRBS[vol] = -1e44;
       for ( nobs=0; nobs<num_obs_angles; nobs++ ) 
       {
         int ind2d=vol+intvolumes*nobs;
         locVolDBS[ind2d] = -1e44;
         locVolDRBS[ind2d] = -1e44;
       }
     }

     if (verbose)
	CCTK_VInfo(CCTK_THORNSTRING,"Calculating total Bremsstrahlung at iteration %d.  Ready for local loops.",cctk_iteration);

     /* Warning: This would break if one patch/level doesn't exist on all cores */
     const int ierri = sumLocalBrems(CCTK_PASS_CTOC, locVolBS, locVolRBS, locVolDBS, locVolDRBS);
     assert(!ierri);
     if (ierri) CCTK_WARN (CCTK_WARN_ABORT, "Could not sum Bremsstrahlung");

     /* Get the summation operator */
     const int sum_op = CCTK_ReductionHandle ("sum");
     if (sum_op < 0)  CCTK_WARN (CCTK_WARN_ABORT, "error");

     /* do all reductions at once to reduce communication overhead if possible */
     /* ReduceArrays is currently only documented inside src/comm/Reduction.c;
      * it is the function underlying ReduceArray. CarpetReduce's reduce.cc
      * reveals that a) num_out_vals is per input array and that num_outvals==1
      * triggers a local reduction whereas b) num_outval==sifzeof(in_array)
      * does an element-by-element reduction which is what we want */
     enum out_buffer_indices { sumBS=0, sumRBS, NUM_BSTYPES };
     CCTK_REAL in_buffer_BS[intvolumes][num_obs_angles+1];
     CCTK_REAL in_buffer_RBS[intvolumes][num_obs_angles+1];
     CCTK_REAL out_buffer[NUM_BSTYPES][intvolumes][num_obs_angles+1];
     for ( vol=0; vol<intvolumes; vol++ )
     {
         in_buffer_BS[vol][0]=locVolBS[vol];
         in_buffer_RBS[vol][0]=locVolRBS[vol];
         out_buffer[sumBS][vol][0]=-1e45;
         out_buffer[sumRBS][vol][0]=-1e45;
         for ( nobs=0; nobs<num_obs_angles; nobs++ ) /* 0 taken by unbeamed */
         {
             int ind2d=vol+intvolumes*nobs;
             in_buffer_BS[vol][nobs+1]=locVolDBS[ind2d];
             in_buffer_RBS[vol][nobs+1]=locVolDRBS[ind2d];
             out_buffer[sumBS][vol][nobs+1]=-1e45;
             out_buffer[sumRBS][vol][nobs+1]=-1e45;
         }
     }
     /* Do reduction */
     ierr = CCTK_ReduceArray (cctkGH,
                              0,                    /* destination processor */
                              sum_op,               /* operation_handle */
                              intvolumes*(num_obs_angles+1), /* num_out_vals */
                              CCTK_VARIABLE_REAL,   /* type_out_vals */
                              out_buffer,           /* out_vals */ 
                              1,                    /* num_dims */
                              NUM_BSTYPES,          /* num_in_fields */
                              CCTK_VARIABLE_REAL,   /* type_in-arrays, last typed argument */
                              intvolumes*(num_obs_angles+1), /* dim[0] (int[num_dims], vararg) */
                              in_buffer_BS,  /* in_arrays, (void *[num_in_fiels], varargs) */
                              in_buffer_RBS );
     assert(!ierr);
     if (ierr) CCTK_WARN (CCTK_WARN_ABORT, "Error in ReduceArrays");

     for (vol=0; vol<intvolumes; vol++)
     {
     	TotBrems[vol] += *SymmetryFactor * out_buffer[sumBS][vol][0];
     	TotRelBrems[vol] += *SymmetryFactor * out_buffer[sumRBS][vol][0];
        for ( nobs=0; nobs<num_obs_angles; nobs++ )
        {
          int ind2d=vol+intvolumes*nobs;
     	  TotDopBrems[ind2d] += out_buffer[sumBS][vol][nobs+1];
     	  TotDopRelBrems[ind2d] += out_buffer[sumRBS][vol][nobs+1];
        }
        if (verbose) {
	   CCTK_VInfo(CCTK_THORNSTRING,"Multiplied level results by %d for Symmetries. ",*SymmetryFactor);
	   CCTK_VInfo(CCTK_THORNSTRING,"Within r=%g, added quantities (dLB,dLRB,dLBth,dLRBth)=(%g,%g,%g,%g)", 
		int_wi_rad[vol], *SymmetryFactor * out_buffer[sumBS][vol][0],  *SymmetryFactor * out_buffer[sumRBS][vol][0], 
		out_buffer[sumBS][vol][1], out_buffer[sumRBS][vol][1] );
        }
     }

  }

}
     
static void bremsstrahlung_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  assert(c >= 0);
  assert(c < intvolumes);
  assert(xc != 0);
  assert(yc != 0);
  assert(zc != 0);

  if (!CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    if (surface_index[c] != -1)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Parameter %s::surface_index[%d] has been set, but centre %d "
                  "is not being set from a spherical surface, so this parameter will be "
                  "ignored.", CCTK_THORNSTRING, c, c);
    }
  }

  if (CCTK_EQUALS(centre_from[c], "parameter"))
  {
    *xc = centre_x[c]; *yc = centre_y[c]; *zc = centre_z[c];
  }
  else if (CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    int sn=surface_index[c];
    if (sn < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d is invalid for volume %d", sn,c);
    } else if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 sn,sphericalsurfaces_nsurfaces);
    }
    if (sf_valid[sn]<=0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "didn't find valid volume centre for sn=%d, vol=%d",sn,c);
    }

    *xc = sf_centroid_x[surface_index[c]];
    *yc = sf_centroid_y[surface_index[c]];
    *zc = sf_centroid_z[surface_index[c]];
  }
  else
  {
    assert(0); // Should never get here
  }

  if (verbose > 0)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "centre %d, at (%g,%g,%g) radius %g", c, *xc, *yc, *zc, int_wi_rad[c]);
  }
}

static int sumLocalBrems(CCTK_ARGUMENTS, CCTK_REAL * locVolBS, 
   CCTK_REAL * locVolRBS, CCTK_REAL * locVolDBS, CCTK_REAL * locVolDRBS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL BremsL,RelBremsL,weightL,rhoL;
  CCTK_REAL gxxL,gyyL,gzzL,gxyL,gxzL,gyzL,detg;
  CCTK_REAL vxL,vyL,vzL,v2,w2L,vcosobs;
  CCTK_REAL xL,yL,zL,r2L;
  CCTK_REAL x0L[intvolumes],y0L[intvolumes],z0L[intvolumes];
  CCTK_REAL dx,dy,dz,dV,dVFact;
  CCTK_REAL intBremsL, intRelBremsL;
  CCTK_REAL Df4[num_obs_angles];
  CCTK_REAL *VolWeights[intvolumes];
  int vol, nobs, my_calc_bs_every;

  for (vol=0; vol<intvolumes; vol++)
  {
    /* Make sure things are initialized */
    locVolBS[vol] = 0;
    locVolRBS[vol] = 0;
    for ( nobs=0; nobs<num_obs_angles; nobs++ ) {
        int ind2d=vol+intvolumes*nobs;
        locVolDBS[ind2d] = 0;
        locVolDRBS[ind2d] = 0;
    }

    /* Get Weight Pointers */
    VolWeights[vol] = NULL;
    VolWeights[vol] = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_var[vol]));
    assert(VolWeights[vol]);
    /* set up local origins for individual integration volumes */
    bremsstrahlung_get_centre(CCTK_PASS_CTOC, vol, &x0L[vol], &y0L[vol], &z0L[vol]);
  }

  /* Summation Limits */
  CCTK_INT   is_symbnd[6], is_physbnd[6], is_ipbnd[6];
  CCTK_INT   nboundaryzones[6], is_internal[6], is_staggered[6], shiftout[6];
  CCTK_INT   imin[3], imax[3], npoints;

  GenericFD_GetBoundaryInfo(cctkGH, cctk_ash, cctk_lsh, cctk_bbox,
                            cctk_nghostzones,
                            imin, imax, is_symbnd, is_physbnd, is_ipbnd);

  if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
    int ierr = GetBoundarySpecification(6, nboundaryzones, is_internal, is_staggered, shiftout);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
  } else {
    /* assume that nboundaryzones = nghostzones
     * leave is_internal, is_staggered, shiftout undefined
     */
    for (int dir = 0; dir<3; dir++)
    {
      for (int face = 0; face < 2; face++)
      {
         CCTK_INT indx = dir*2 + face;
         nboundaryzones[indx] = cctk_nghostzones[dir];
      }
    }
  }

  for (int dir = 0; dir<3; dir++)
  {
    for (int face = 0; face < 2; face++)
    {
       /* Include physical boundaries, not sym bdies or ghostzones */
       CCTK_INT indx = dir*2 + face;
       if (is_physbnd[indx]) {
           npoints=0;
       } else if (is_symbnd[indx]) {
           npoints=nboundaryzones[indx];
       } else {
           npoints=cctk_nghostzones[dir];
       }

       switch(face)
       {
       case 0: /*lower*/
          imin[dir]=npoints;
          break;
       case 1: /*upper*/
          imax[dir]=cctk_lsh[dir] - npoints;
          break;
       default:
          CCTK_WARN(0,"internal error");
       }

    }
  }

  if (verbose) { 
	CCTK_VInfo(CCTK_THORNSTRING,"Extent: (%d,%d,%d) to (%d,%d,%d)",imin[0],imin[1],imin[2],imax[0],imax[1],imax[2]);
  }


  /* Loop over everything except ghost zones */ 
  for(int k = imin[2] ; k < imax[2] ; k++)
  {
    for(int j = imin[1] ; j < imax[1] ; j++)
    {
      for(int i = imin[0] ; i < imax[0] ; i++)
      {

        /* Initiate variables */
        BremsL = INITVALUE;
        RelBremsL = INITVALUE;
        weightL = INITVALUE;
	gxxL = INITVALUE;
	gxyL = INITVALUE;
	gxzL = INITVALUE;
	gyyL = INITVALUE;
	gyzL = INITVALUE;
	gzzL = INITVALUE;

        dx = INITVALUE;
        dy = INITVALUE;
        dz = INITVALUE;
        dV = INITVALUE;
	xL = INITVALUE;
	yL = INITVALUE;
	zL = INITVALUE;
	r2L = INITVALUE;

	vxL = INITVALUE;
	vyL = INITVALUE;
	vzL = INITVALUE;
	v2 = INITVALUE;
	w2L = INITVALUE;

        intBremsL = INITVALUE;
        intRelBremsL = INITVALUE;

	/* Information from the grid */
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);

	rhoL = rho[idx];
	BremsL = Brems[idx];
	RelBremsL = RelBrems[idx];

	xL = x[idx];
	yL = y[idx];
	zL = z[idx];

	gxxL = gxx[idx];
	gxyL = gxy[idx];
	gxzL = gxz[idx];
	gyyL = gyy[idx];
	gyzL = gyz[idx];
	gzzL = gzz[idx];

#ifdef HAVE_MAYA_WHISKY
	vxL = velx[idx];
	vyL = vely[idx];
	vzL = velz[idx];
#else
        int xidx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, 0);
        int yidx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, 1);
        int zidx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k, 2);
        vxL = vel[xidx]; 
        vyL = vel[yidx];
        vzL = vel[zidx];
#endif
	if (debug)
	{
	   vyL = 0.3;
	   vzL = 0.4;
	}

	v2 = gxxL*SQR(vxL) + gyyL*SQR(vyL) + gzzL*SQR(vzL) 
	     + 2*gxyL*vxL*vyL + 2*gxzL*vxL*vzL * 2*gyzL*vyL*vzL;
	w2L = 1. / (1. - v2);

	if (debug)
	   w2L = 4./3.;

	/* Doppler factors, observer at -\hat{y} and varying angles above plane */
        CCTK_REAL tanth, costh, obsmag;
        for ( nobs=0; nobs<num_obs_angles; nobs++ ) {
            Df4[nobs] = 0;
            tanth=tan( obs_angles[nobs] );
            obsmag=sqrt( 1. + SQR(tanth) );
            for ( int term=0; term<n_DTerms; term++ ) {
	        CCTK_REAL vxLL = Sym_Signs[3*term]*vxL;
	        CCTK_REAL vyLL = Sym_Signs[3*term+1]*vyL;
                CCTK_REAL vzLL = Sym_Signs[3*term+2]*vzL;
		if ( abs(tanth) <= 1. ) {
                   vcosobs = (tanth*vzLL - vyLL)/obsmag; /* theta rotates observer in yz-plane. Avoid the poles! */
		} else {
		   tanth=tan( PI/2. - obs_angles[nobs] ); /* Switch theta to measure from pole */
                   obsmag=sqrt( 1. + SQR(tanth) );
                   vcosobs = (vzLL - tanth*vyLL)/obsmag;
		}
                Df4[nobs] += SQR(1./(w2L * SQR(1. - vcosobs)));
                if ( idx%1000==0 && verbose ) {
                   CCTK_VInfo(CCTK_THORNSTRING,"For term %d, Signs=(%d,%d,%d), vLL=(%g,%g,%g), vcosobs=%g, Df4+=%g.",
		   term,Sym_Signs[3*term],Sym_Signs[3*term+1],Sym_Signs[3*term+2],vxLL,vyLL,vzLL,vcosobs, 
		   SQR(1./(w2L * SQR(1. - vcosobs))));
	        }
            }
        }
        if ( idx%1000==0 && verbose ) {
           for ( nobs=0; nobs<num_obs_angles; nobs++ ) {
            CCTK_VInfo(CCTK_THORNSTRING," ... For observer %d, Df4=%g.",nobs,Df4[nobs]);
           }
	}

	/* Now for integration .... */
        /* Get volume element       */
        dx = CCTK_DELTA_SPACE(0);
        dy = CCTK_DELTA_SPACE(1);
        dz = CCTK_DELTA_SPACE(2);

	/* detg */
	detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
	detg = sqrt(detg);
	if (do_cgs_units) {
      		dVFact=3.24e15 * Msys_SolarMasses*SQR(Msys_SolarMasses);
  	} else {
      		dVFact=1;
  	}
        dV = dVFact * detg * (dx*dy*dz);

	/* Debug */
	if (debug) {
		BremsL = 1;
		RelBremsL = 0.01;
		dV = dx*dy*dz; /* Override phys_int */
	}

	if ( rhoL <= rho_min ) 
	{
		continue;
		intBremsL = 0;
		intRelBremsL = 0;
	} else {
		intBremsL = BremsL*dV;
		intRelBremsL = RelBremsL*dV;
	}


	/* Calculate integrand */ 
	for (vol=0; vol<intvolumes; vol++)
	{
           my_calc_bs_every = calc_bs_every_vol[vol] >= 0 ? calc_bs_every_vol[vol] : calc_bs_every;
           assert(my_calc_bs_every >= 0);
           if ( (my_calc_bs_every == 0) || (cctk_iteration % my_calc_bs_every != 0) )
               continue;

	   r2L = SQR(xL-x0L[vol]) + SQR(yL-y0L[vol]) + SQR(zL-z0L[vol]) ;
	   if ( r2L < SQR(int_wi_rad[vol]) )
	   {
              weightL = VolWeights[vol][idx];
              locVolBS[vol]   += weightL*intBremsL;
              locVolRBS[vol]  += weightL*intRelBremsL;
              for ( nobs=0; nobs<num_obs_angles; nobs++ ) {
                  int ind2d=vol+intvolumes*nobs;
                  locVolDBS[ind2d]  += weightL*Df4[nobs]*intBremsL;
                  locVolDRBS[ind2d] += weightL*Df4[nobs]*intRelBremsL;
              }
	   }
	}
        
      }
    }
  }

  if (verbose){
     CCTK_INFO("After local integration, local totals are:");
     for (vol=0; vol<intvolumes; vol++)
     {
       CCTK_VInfo(CCTK_THORNSTRING,"\t Within r=%g: Brems = %.17g, RelBrems = %.17g, DopBrems = %.17g, DopRelBrems = %.17g", 
	          int_wi_rad[vol], locVolBS[vol], locVolRBS[vol], locVolDBS[vol], locVolDRBS[vol]);
     }
  }

  return 0;

}

/* IO */
void Brems_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const char *fmode = (file_created>0) ? "a" : "w";  // file mode: append if already written
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  char addformat[50];
  int vol, nobs, do_calc_bs, my_calc_bs_every;
  FILE *file;

  const CCTK_INT myproc= CCTK_MyProc(cctkGH);

  do_calc_bs = 0;
  for (vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_bs_every = calc_bs_every_vol[vol] >= 0 ? calc_bs_every_vol[vol] : calc_bs_every;
    assert(my_calc_bs_every >= 0);
    if ( (my_calc_bs_every > 0) && (cctk_iteration % my_calc_bs_every == 0) )
    {
        do_calc_bs = 1;
        break;
    }
  }
  if ( (myproc == 0) && do_calc_bs ) 
  {

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "writing output");
    }

    // filename
    sprintf(varname, "BSLuminosities");

    filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                                strlen (file_extension) +2);
    assert(filename);
    sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

    // open file
    file = fopen (filename, fmode);
    if (!file) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "write_Brems: Could not open integrated scalar output file '%s'",
                  filename);
    }

    // write header on startup
    if (file_created<=0) {
      fprintf(file,"# Integrated Bremsstrahlung luminosities.\n");
      fprintf(file,"# gnuplot column index:\n");
      fprintf(file,"# 1:it 2:t \n");
      int LperV=2*(num_obs_angles+1);
      for (vol=0;vol<intvolumes;vol++)
      {
	fprintf(file,"#\t( r=%g ) %d:L_BS %d:L_RBS ",int_wi_rad[vol],LperV*vol+3,LperV*vol+4);
        for ( nobs=0; nobs<num_obs_angles; nobs++ )
        {
	    fprintf(file," %d:L_BS(th=%3.4g) %d:L_RBS(th=%3.4g)",LperV*vol+2*nobs+5,obs_angles[nobs],LperV*vol+2*nobs+6,obs_angles[nobs]);
        }
        fprintf(file,"\n");
      }
    }

    // format string
    sprintf (format_str_real,"%%d\t%%%s",out_format);
    sprintf (addformat,"\t%%%s\t%%%s",out_format,out_format);
    // write data
    fprintf(file, format_str_real, cctk_iteration, cctk_time );
    for (vol=0;vol<intvolumes;vol++)
    {
	fprintf(file,addformat,TotBrems[vol],TotRelBrems[vol]);
        for (nobs=0; nobs<num_obs_angles; nobs++ )
        {
          int ind2d=vol+intvolumes*nobs;
	  fprintf(file,addformat,TotDopBrems[ind2d],TotDopRelBrems[ind2d]);
        }
    }
    fprintf(file,"\n");

    fclose(file);
    free(filename);

    if (file_created==0) {
      file_created=1;
    }
  } 
}
