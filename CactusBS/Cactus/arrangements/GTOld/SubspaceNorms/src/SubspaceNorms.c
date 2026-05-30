#define INITVALUE (42)

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define KRANC_C
#include "GenericFD.h"
#include "Symmetry.h"

#define SQR(x) ((x) * (x))
#define MAX_INTVOLUMES 100
#define MAX_NUMBER_NORMS 20

//*******************************//
//****** External Routines ******//
//*******************************//
void SubspaceNorms_Scalars(CCTK_ARGUMENTS);
void SubspaceNorms_Init(CCTK_ARGUMENTS);
void ZeroNorms(CCTK_ARGUMENTS);
void SumSubspaceNorms(CCTK_ARGUMENTS);
void SubspaceNorms_Output(CCTK_ARGUMENTS);

//*******************************//
//****** Internal Routines ******//
//*******************************//
static int sumNormsLocally(CCTK_ARGUMENTS, CCTK_REAL * localNorms, CCTK_REAL * localCount);
static void subspacenorms_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc, CCTK_REAL *radc);
static void fill_var_idx_array(int idx, const char *optstring, void *callback_arg);
static void fill_grp_idx_array(int idx, const char *optstring, void *callback_arg);

//******************************//
//****** Global Variables ******//
//******************************//
static CCTK_INT file_created[MAX_INTVOLUMES];


//***********************//
//****** Functions ******//
//***********************//

void SubspaceNorms_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT   is_symbnd[6];

  *SymmetryFactor = INITVALUE;

  for ( int vol=0; vol<intvolumes; vol++) {
      file_created[vol]=0;
      for ( int nscalar=0; nscalar<num_scalars; nscalar++ ) {
          ScalarNorms[ vol*num_scalars + nscalar ] = -1e234;
      }
      for ( int nvector=0; nvector<num_vectors; nvector++ ) {
          VectorNorms[ vol*num_vectors + nvector ] = -1e234;
      }
      for ( int ntensor=0; ntensor<num_tensors; ntensor++ ) {
          TensorNorms[ vol*num_tensors + ntensor ] = -1e234;
      }
      for ( int norm_in_var=0; norm_in_var<(num_scalars+num_vectors+num_tensors); norm_in_var++ ) {
          NormCounter[ vol*(num_scalars+num_vectors+num_tensors) + norm_in_var ] = -1e234;
      }
  }

  /* Get symmetry info */
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
  /*********************/

  /** Initialize index arrays, then fill **/
  for ( int idx=0; idx<num_scalars; idx++ ) {
      varinfo_scalars[idx] = -1;
  }
  for ( int idx=0; idx<num_vectors; idx++ ) {
      groupinfo_vectors[idx] = -1;
  }
  for ( int idx=0; idx<num_tensors; idx++ ) {
      groupinfo_tensors[idx] = -1;
  }

  if ( !CCTK_EQUALS(scalar_variables, "")) {
     ierr = CCTK_TraverseString( scalar_variables, fill_var_idx_array, varinfo_scalars, CCTK_VAR );
     assert(ierr>0);
     if ( verbose ) {
        for ( int nscalar=0; nscalar<num_scalars; nscalar++ ) {
	    CCTK_VInfo(CCTK_THORNSTRING,"Asked for norms of variable %s, treated as a scalar.", CCTK_VarName(varinfo_scalars[nscalar]));
	    
        } 
     }
  }

  if ( !CCTK_EQUALS(vector_groups, "")) {
     ierr = CCTK_TraverseString( vector_groups, fill_grp_idx_array, groupinfo_vectors, CCTK_GROUP );
     assert(ierr>0);
     for ( int nvector=0; nvector<num_vectors; nvector++ ) {
        if ( verbose )
	    CCTK_VInfo(CCTK_THORNSTRING,"Asked for norms of group %s, treated as a vector.", CCTK_GroupName(groupinfo_vectors[nvector]));
        if ( CCTK_NumVarsInGroupI( groupinfo_vectors[nvector] ) != 3 )
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Asked for vector treatment of a non-vector %s (at least, num_vars(group) =/= 3)   ", CCTK_GroupName(groupinfo_vectors[nvector]) );
     }
  }

  if ( !CCTK_EQUALS(tensor_groups, "")) {
     ierr = CCTK_TraverseString( tensor_groups, fill_grp_idx_array, groupinfo_tensors, CCTK_GROUP );
     assert(ierr>0);
     for ( int ntensor=0; ntensor<num_tensors; ntensor++ ) {
         if ( verbose )
	    CCTK_VInfo(CCTK_THORNSTRING,"Asked for norms of group %s, treated as a tensor.", CCTK_GroupName(groupinfo_tensors[ntensor]));
         if ( CCTK_NumVarsInGroupI( groupinfo_tensors[ntensor] ) != 6 )
            CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Asked for tensor treatment of a non-tensor %s (at least, num_vars(group) =/= 6)   ", CCTK_GroupName(groupinfo_tensors[ntensor]) ); 
     }
  }

}

/* Callback routine to put the cactus var index into a given array. */
static void fill_var_idx_array(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg);

  CCTK_INT *idx_array = (CCTK_INT * ) callback_arg;

  /* find the first free slot and store the new index in it */
  for(int i = 0 ; i < MAX_NUMBER_NORMS ; i++)
  {
    if(idx_array[i] == -1) {
      idx_array[i] = idx;
      break;
    }
  }
}

/* Callback routine to put the cactus var index into a given array. */
static void fill_grp_idx_array(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg);

  CCTK_INT *idx_array = (CCTK_INT * ) callback_arg;
  // Traverse string gives the variable indicies of all variables in group //
  int grpidx = CCTK_GroupIndexFromVarI(idx);
  if ( idx != CCTK_FirstVarIndexI( grpidx ) )
     return;

  /* find the first free slot and store the new index in it */
  for(int i = 0 ; i < MAX_NUMBER_NORMS ; i++)
  {
    if(idx_array[i] == -1) {
      idx_array[i] = grpidx;
      break;
    }
  }
}

void ZeroNorms(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for ( int vol=0; vol<intvolumes; vol++) {

      for ( int nscalar=0; nscalar<num_scalars; nscalar++ )
          ScalarNorms[ vol*num_scalars + nscalar ] = 0;

      for ( int nvector=0; nvector<num_vectors; nvector++ )
          VectorNorms[ vol*num_vectors + nvector ] = 0;

      for ( int ntensor=0; ntensor<num_tensors; ntensor++ )
          TensorNorms[ vol*num_tensors + ntensor ] = 0;

      for ( int norm_in_vol=0; norm_in_vol<(num_scalars+num_vectors+num_tensors); norm_in_vol++ ) 
          NormCounter[ vol*(num_scalars+num_vectors+num_tensors) + norm_in_vol ] = 0;
  }
}


void SumSubspaceNorms(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT do_calc_norms, my_calc_norms_every;
  CCTK_INT norms_per_vol = num_scalars + num_vectors + num_tensors;

  CCTK_REAL localNorms[intvolumes*norms_per_vol];
  CCTK_REAL localCounts[intvolumes*norms_per_vol];
  CCTK_REAL sumNorms[intvolumes*norms_per_vol];
  CCTK_REAL sumCounts[intvolumes*norms_per_vol];

 
  do_calc_norms = 0;
  for ( int vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_norms_every = calc_norms_every_vol[vol] >= 0 ? calc_norms_every_vol[vol] : calc_norms_every;
    assert(my_calc_norms_every >= 0);
    if ( (my_calc_norms_every > 0) && (cctk_iteration % my_calc_norms_every == 0) )
    {
	/* Go through loop if any single volume is to be calculated */
        do_calc_norms = 1;
        break;
    }
  }

  if (do_calc_norms) 
  {

     if (verbose)
	CCTK_VInfo(CCTK_THORNSTRING,"Calculating subregion norms at iteration %d",cctk_iteration);

     /* Initialize local versions for local Norms */
     for ( int vol=0; vol<intvolumes; vol++)
     {
         for ( int norm_in_vol=0; norm_in_vol<norms_per_vol; norm_in_vol++ ) {
             localNorms [ vol*norms_per_vol + norm_in_vol ] = -1e44;
             localCounts[ vol*norms_per_vol + norm_in_vol ] = -1e44;
             sumNorms [ vol*norms_per_vol + norm_in_vol ] = 0;
             sumCounts[ vol*norms_per_vol + norm_in_vol ] = 0;
         }
     }

     if (verbose) CCTK_INFO("Calling local loops for processor-local integrations.");
     const int ierrL = sumNormsLocally(CCTK_PASS_CTOC, localNorms, localCounts);
     assert(!ierrL);
     if (ierrL) CCTK_WARN (CCTK_WARN_ABORT, "Error in performing processor-local integrations.");

     const CCTK_INT myproc = CCTK_MyProc(cctkGH);
     if (verbose) { /* Dump all local data to standard out */
        CCTK_VInfo( CCTK_THORNSTRING, "Local integrations on processor %d, iteration %d:", myproc, cctk_iteration );
        for ( int vol=0; vol<intvolumes; vol++ ) 
            for ( int norm=0; norm<norms_per_vol; norm++ )
                CCTK_VInfo( CCTK_THORNSTRING, "    (%d,%d,%d): Norm=%g, Count=%g", myproc, vol, 
                      norm, localNorms[vol*norms_per_vol+norm], localCounts[vol*norms_per_vol+norm] );
     }

     /* Get the summation operator */
     const int sum_op = CCTK_ReductionHandle ("sum");
     if (sum_op < 0)  CCTK_WARN (CCTK_WARN_ABORT, "error");

     /* Currently doing reductions in 2 steps. */
     int ierrN = CCTK_ReduceArray (cctkGH, 0, sum_op, intvolumes*norms_per_vol, CCTK_VARIABLE_REAL, sumNorms, 
                                   1,1, CCTK_VARIABLE_REAL, intvolumes*norms_per_vol, localNorms);
     assert(!ierrN);
     if (ierrN) CCTK_WARN (CCTK_WARN_ABORT, "Error in reducing localNorms across processors.");

     int ierrC = CCTK_ReduceArray (cctkGH, 0, sum_op, intvolumes*norms_per_vol, CCTK_VARIABLE_REAL, sumCounts, 
                                   1,1, CCTK_VARIABLE_REAL, intvolumes*norms_per_vol, localCounts);
     assert(!ierrC);
     if (ierrC) CCTK_WARN (CCTK_WARN_ABORT, "Error in reducing Counters across processors.");

     if ( myproc==0 ) {

        if (verbose) {
           CCTK_VInfo( CCTK_THORNSTRING, "Merged integrations on processor %d, iteration %d:", myproc, cctk_iteration );
           for ( int vol=0; vol<intvolumes; vol++ ) 
               for ( int norm=0; norm<norms_per_vol; norm++ )
                   CCTK_VInfo( CCTK_THORNSTRING, "    (%d,%d): Norm=%g, Count=%g", vol, norm, 
                         sumNorms[vol*norms_per_vol+norm], sumCounts[vol*norms_per_vol+norm] );
        }
   
        /* Write to grid arrays, applying symmetry factor       */
        /*   (we add because this entire function is called at  */
        /*          each refinement level separately)           */
        for ( int vol=0; vol<intvolumes; vol++) {
            for ( int nscalar=0; nscalar<num_scalars; nscalar++ ) {
                ScalarNorms[ vol*num_scalars + nscalar ] += *SymmetryFactor*sumNorms[ vol*norms_per_vol + nscalar ];
            }
            for ( int nvector=0; nvector<num_vectors; nvector++ ) {
                VectorNorms[ vol*num_vectors + nvector ] += *SymmetryFactor*sumNorms[ vol*norms_per_vol + num_scalars + nvector ];
            }
            for ( int ntensor=0; ntensor<num_tensors; ntensor++ ) {
                TensorNorms[ vol*num_tensors + ntensor ] += *SymmetryFactor*sumNorms[ vol*norms_per_vol + num_scalars + num_vectors + ntensor ];
            }
            for ( int norm_in_vol=0; norm_in_vol<norms_per_vol; norm_in_vol++ ) {
                NormCounter[ vol*norms_per_vol + norm_in_vol ] += (*SymmetryFactor)*sumCounts[ vol*norms_per_vol + norm_in_vol ];
            }
        }
        
        if (verbose) 
   	   CCTK_VInfo(CCTK_THORNSTRING,"Multiplied level results by %d for Symmetries: ",*SymmetryFactor);

     }

  }

}
     
static void subspacenorms_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc, CCTK_REAL *radc)
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
    *xc = centre_x[c]; *yc = centre_y[c]; *zc = centre_z[c]; *radc = int_wi_rad[c];
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
	*radc = sf_mean_radius[surface_index[c]];
  }
  else
  {
    assert(0); // Should never get here
  }

  if (verbose > 0)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "centre %d, at (%g,%g,%g) radius %g", c, *xc, *yc, *zc, *radc);
  }
}

static int sumNormsLocally(CCTK_ARGUMENTS, CCTK_REAL * localNorms, CCTK_REAL * localCount)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL x0L[intvolumes],y0L[intvolumes],z0L[intvolumes],rad0L[intvolumes];

  /* Pointers to actual Gridfunction data */
  CCTK_REAL *ScalarVariable[num_scalars];
  CCTK_REAL *ScalarFloorVar[num_scalars];
  CCTK_REAL *VectorVariable_x[num_vectors];
  CCTK_REAL *VectorVariable_y[num_vectors];
  CCTK_REAL *VectorVariable_z[num_vectors];
  CCTK_REAL *TensorVariable_xx[num_tensors];
  CCTK_REAL *TensorVariable_xy[num_tensors];
  CCTK_REAL *TensorVariable_xz[num_tensors];
  CCTK_REAL *TensorVariable_yy[num_tensors];
  CCTK_REAL *TensorVariable_yz[num_tensors];
  CCTK_REAL *TensorVariable_zz[num_tensors];
  CCTK_REAL *VolWeights[intvolumes];
  CCTK_REAL *ExtraWeightVar[num_scalars];

  /* Store values at a given point */
  CCTK_REAL scalarVal[num_scalars];
  CCTK_REAL ScalarFloorVal[num_scalars];
  CCTK_REAL ExtraWeightVal[num_scalars];
  CCTK_REAL vectorVal_x[num_vectors];
  CCTK_REAL vectorVal_y[num_vectors];
  CCTK_REAL vectorVal_z[num_vectors];
  CCTK_REAL tensorVal_xx[num_tensors];
  CCTK_REAL tensorVal_xy[num_tensors];
  CCTK_REAL tensorVal_xz[num_tensors];
  CCTK_REAL tensorVal_yy[num_tensors];
  CCTK_REAL tensorVal_yz[num_tensors];
  CCTK_REAL tensorVal_zz[num_tensors];

  CCTK_INT norms_per_vol = num_scalars+num_vectors+num_tensors;

  /* Make sure things are initialized */
  for ( int vol=0; vol<intvolumes; vol++)
  {
    /* Make sure local variables are initialized */
    for ( int norm_in_vol=0; norm_in_vol<norms_per_vol; norm_in_vol++ ) {
        localNorms[ vol*norms_per_vol+norm_in_vol ] = 0;
        localCount[ vol*norms_per_vol+norm_in_vol ] = 0;
    }

    /* Get Weight Pointers */
    VolWeights[vol] = NULL;
    VolWeights[vol] = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_var[vol]));
    assert(VolWeights[vol]);

    /* set up local origins for individual integration volumes */
    subspacenorms_get_centre(CCTK_PASS_CTOC, vol, &x0L[vol], &y0L[vol], &z0L[vol], &rad0L[vol]);

  }

  /* Get Pointers to Data */
  for ( int nscalar=0; nscalar<num_scalars; nscalar++ ) {

      ScalarVariable[nscalar]=NULL;
      ScalarVariable[nscalar]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0,varinfo_scalars[nscalar]));
      assert(ScalarVariable[nscalar]);
      if (verbose) {
         CCTK_VInfo(CCTK_THORNSTRING, "Pointing ScalarVariable[%d] to %s", nscalar,CCTK_VarName( varinfo_scalars[nscalar] ));

      }
      /* Get Scalar Floor Pointers */
      ScalarFloorVar[nscalar] = ScalarVariable[nscalar];
      CCTK_INT ScalarFloorIdx = varinfo_scalars[nscalar];
      if ( scalar_floor[nscalar] ) {
         if ( strcmp(scalar_floor_var[nscalar]," ") != 0 ) {
            ScalarFloorVar[nscalar] = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,scalar_floor_var[nscalar]));
            ScalarFloorIdx = CCTK_VarIndex( scalar_floor_var[nscalar] );
            assert(ScalarFloorVar[nscalar]);
         }
         if (verbose)
            CCTK_VInfo(CCTK_THORNSTRING, "Pointing ScalarFloorVar[%d] to %s", nscalar, CCTK_VarName(ScalarFloorIdx)); 
      }

      /* Scalar-associated extra density weight factor */
      if ( density_weighted[nscalar] > 0 ) {
         ExtraWeightVar[nscalar] = NULL;
         ExtraWeightVar[nscalar] = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,"HydroBase::rho"));
         assert(ExtraWeightVar[nscalar]);
      }

  }

  for ( int nvector=0; nvector<num_vectors; nvector++ ) {
      VectorVariable_x[nvector]=0;
      VectorVariable_y[nvector]=0;
      VectorVariable_z[nvector]=0;

      int first_varIdx = CCTK_FirstVarIndexI( groupinfo_vectors[nvector] );
      VectorVariable_x[nvector]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx));
      VectorVariable_y[nvector]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+1));
      VectorVariable_z[nvector]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+2));

      assert(VectorVariable_x[nvector]);
      assert(VectorVariable_y[nvector]);
      assert(VectorVariable_z[nvector]);
  }

  for ( int ntensor=0; ntensor<num_tensors; ntensor++ ) {
      TensorVariable_xx[ntensor]=0;
      TensorVariable_xy[ntensor]=0;
      TensorVariable_xz[ntensor]=0;
      TensorVariable_yy[ntensor]=0;
      TensorVariable_yz[ntensor]=0;
      TensorVariable_zz[ntensor]=0;

      int first_varIdx = CCTK_FirstVarIndexI( groupinfo_tensors[ntensor] );
      TensorVariable_xx[ntensor]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx));
      TensorVariable_xy[ntensor]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+1));
      TensorVariable_xz[ntensor]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+2));
      TensorVariable_yy[ntensor]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+3));
      TensorVariable_yz[ntensor]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+4));
      TensorVariable_zz[ntensor]=(CCTK_REAL *)(CCTK_VarDataPtrI(cctkGH,0, first_varIdx+5));

      assert(TensorVariable_xx[ntensor]);
      assert(TensorVariable_xy[ntensor]);
      assert(TensorVariable_xz[ntensor]);
      assert(TensorVariable_yy[ntensor]);
      assert(TensorVariable_yz[ntensor]);
      assert(TensorVariable_zz[ntensor]);
  }

  /* Summation limits */
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
//  #pragma omp parallel for reduction ( +:localNorms,localCounts ) 
  for(int k = imin[2] ; k < imax[2] ; k++)
  {
    for(int j = imin[1] ; j < imax[1] ; j++)
    {
      for(int i = imin[0] ; i < imax[0] ; i++)
      {

        /* Declare here due to openmp */
        CCTK_REAL gxxL,gyyL,gzzL,gxyL,gxzL,gyzL,detg;
        CCTK_REAL gInv_xx,gInv_yy,gInv_zz,gInv_xy,gInv_xz,gInv_yz;
        CCTK_REAL dx,dy,dz,dV,weightL;
        CCTK_REAL xL,yL,zL,r2L;
        CCTK_REAL dNorm[norms_per_vol];

        /* Initiate variables */
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

	for( int nscalar=0 ; nscalar<num_scalars ; nscalar++ ) {
	   scalarVal[nscalar] = INITVALUE;
           ScalarFloorVal[nscalar] = INITVALUE;
	}
	for( int nvector=0 ; nvector<num_vectors ; nvector++ ) {
	   vectorVal_x[nvector] = INITVALUE;
	   vectorVal_y[nvector] = INITVALUE;
	   vectorVal_z[nvector] = INITVALUE;
	}
	for( int ntensor=0 ; ntensor<num_tensors ; ntensor++ ) {
	   tensorVal_xx[ntensor] = INITVALUE;
	   tensorVal_xy[ntensor] = INITVALUE;
	   tensorVal_xz[ntensor] = INITVALUE;
	   tensorVal_yy[ntensor] = INITVALUE;
	   tensorVal_yz[ntensor] = INITVALUE;
	   tensorVal_zz[ntensor] = INITVALUE;
	}

	/* Information from the grid */
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);

	xL = x[idx];
	yL = y[idx];
	zL = z[idx];

	gxxL = gxx[idx];
	gxyL = gxy[idx];
	gxzL = gxz[idx];
	gyyL = gyy[idx];
	gyzL = gyz[idx];
	gzzL = gzz[idx];

        /* Get volume element */
        dx = CCTK_DELTA_SPACE(0);
        dy = CCTK_DELTA_SPACE(1);
        dz = CCTK_DELTA_SPACE(2);
        dV = dx*dy*dz;

	for( int nscalar=0 ; nscalar<num_scalars ; nscalar++ ) {
	   scalarVal[nscalar] = ScalarVariable[nscalar][idx];
	   ScalarFloorVal[nscalar] = ScalarFloorVar[nscalar][idx];
           ExtraWeightVal[nscalar] = ( density_weighted[nscalar] ? ExtraWeightVar[nscalar][idx] : 1 );
	}
	for( int nvector=0 ; nvector<num_vectors ; nvector++ ) {
	   vectorVal_x[nvector] = VectorVariable_x[nvector][idx];
	   vectorVal_y[nvector] = VectorVariable_y[nvector][idx];
	   vectorVal_z[nvector] = VectorVariable_z[nvector][idx];
	}
	for( int ntensor=0 ; ntensor<num_tensors ; ntensor++ ) {
	   tensorVal_xx[ntensor] = TensorVariable_xx[ntensor][idx];
	   tensorVal_xy[ntensor] = TensorVariable_xy[ntensor][idx];
	   tensorVal_xz[ntensor] = TensorVariable_xz[ntensor][idx];
	   tensorVal_yy[ntensor] = TensorVariable_yy[ntensor][idx];
	   tensorVal_yz[ntensor] = TensorVariable_yz[ntensor][idx];
	   tensorVal_zz[ntensor] = TensorVariable_zz[ntensor][idx];
	}


	/* Calculate local contribution to norms */
	detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
        gInv_xx = (-SQR(gyzL) + gyyL*gzzL)/detg;
	gInv_xy = (gxzL*gyzL - gxyL*gzzL)/detg;
        gInv_xz = (-gxzL*gyyL + gxyL*gyzL)/detg;
        gInv_yy = (-SQR(gxzL) + gxxL*gzzL)/detg;
        gInv_yz = (gxyL*gxzL - gxxL*gyzL)/detg;
        gInv_zz = (-SQR(gxyL) + gxxL*gyyL)/detg;

	for( int nscalar=0 ; nscalar<num_scalars ; nscalar++ ) {
           if ( scalar_floor[nscalar] > 0 && (ScalarFloorVal[nscalar] < scalar_floor_val[nscalar]) ) {
              dNorm[ nscalar ] = 0;
           } else { 
              dNorm[ nscalar ] = pow(scalarVal[nscalar],norm_exponent[nscalar]);
           } 
	   if ( undensitize_scalar[nscalar] ) {
              dNorm[ nscalar ] = dNorm[nscalar]/sqrt(detg);
	   }
	}

	for( int nvector=0 ; nvector<num_vectors ; nvector++ )
	   dNorm[ num_scalars + nvector ] = ( gxxL*SQR( vectorVal_x[nvector] ) + gyyL*SQR( vectorVal_y[nvector] )
	       + gzzL*SQR( vectorVal_z[nvector] ) + 2*gxyL*vectorVal_x[nvector]*vectorVal_y[nvector]
	       + 2*gxzL*vectorVal_x[nvector]*vectorVal_z[nvector]
	       + 2*gyzL*vectorVal_y[nvector]*vectorVal_z[nvector] );

	for( int ntensor=0 ; ntensor<num_tensors ; ntensor++ ) {

	   CCTK_REAL tInv_xx, tInv_xy, tInv_xz, tInv_yy, tInv_yz, tInv_zz;
           tInv_xx = tensorVal_xx[ntensor]*SQR( gInv_xx ) + tensorVal_yy[ntensor]*SQR( gInv_xy )+ tensorVal_zz[ntensor]*SQR( gInv_xz )
              + 2*tensorVal_xy[ntensor]*gInv_xx*gInv_xy + 2*tensorVal_xz[ntensor]*gInv_xx*gInv_xz
              + 2*tensorVal_yz[ntensor]*gInv_xy*gInv_xz;
           tInv_yy = tensorVal_xx[ntensor]*SQR( gInv_xy ) + tensorVal_yy[ntensor]*SQR( gInv_yy )+ tensorVal_zz[ntensor]*SQR( gInv_yz )
              + 2*tensorVal_xy[ntensor]*gInv_xy*gInv_yy + 2*tensorVal_xz[ntensor]*gInv_xy*gInv_yz
              + 2*tensorVal_yz[ntensor]*gInv_yy*gInv_yz;
           tInv_zz = tensorVal_xx[ntensor]*SQR( gInv_xz ) + tensorVal_yy[ntensor]*SQR( gInv_yz )+ tensorVal_zz[ntensor]*SQR( gInv_zz )
              + 2*tensorVal_xy[ntensor]*gInv_xz*gInv_yz + 2*tensorVal_xz[ntensor]*gInv_xz*gInv_zz
              + 2*tensorVal_yz[ntensor]*gInv_yz*gInv_zz;
           tInv_xy = tensorVal_xx[ntensor]*gInv_xx*gInv_xy + tensorVal_yy[ntensor]*gInv_xy*gInv_yy + tensorVal_zz[ntensor]*gInv_xz*gInv_yz
              + tensorVal_xy[ntensor]*( gInv_xx*gInv_yy + SQR(gInv_xy) )
              + tensorVal_xz[ntensor]*( gInv_xx*gInv_yz + gInv_xy*gInv_xz )
              + tensorVal_yz[ntensor]*( gInv_xx*gInv_yz + gInv_xy*gInv_xz );
           tInv_xz = tensorVal_xx[ntensor]*gInv_xx*gInv_xz + tensorVal_yy[ntensor]*gInv_xy*gInv_yz + tensorVal_zz[ntensor]*gInv_xz*gInv_zz
              + tensorVal_xy[ntensor]*( gInv_xx*gInv_yz + gInv_xz*gInv_xy )
              + tensorVal_xz[ntensor]*( gInv_xx*gInv_zz + SQR(gInv_xz) )
              + tensorVal_yz[ntensor]*( gInv_xy*gInv_zz + gInv_yz*gInv_xz );
           tInv_yz = tensorVal_xx[ntensor]*gInv_xy*gInv_xz + tensorVal_yy[ntensor]*gInv_yy*gInv_yz + tensorVal_zz[ntensor]*gInv_yz*gInv_zz
              + tensorVal_xy[ntensor]*( gInv_xy*gInv_yz + gInv_xz*gInv_yy )
              + tensorVal_xz[ntensor]*( gInv_xy*gInv_zz + gInv_xz*gInv_yz )
              + tensorVal_yz[ntensor]*( gInv_yy*gInv_zz + SQR(gInv_yz) );

           int tidx = num_scalars+num_vectors + ntensor; 
	   dNorm[tidx] = tInv_xx*tensorVal_xx[ntensor] + 2*tInv_xy*tensorVal_xy[ntensor]
               + 2*tInv_xz*tensorVal_xz[ntensor] + tInv_yy*tensorVal_yy[ntensor]
	       + 2*tInv_yz*tensorVal_yz[ntensor] + tInv_zz*tensorVal_zz[ntensor];
           if ( debug && (idx%100==0) ) {
              CCTK_VInfo(CCTK_THORNSTRING,"If the tensor is the metric, g^ab g_ab should be 3.  It is %g", dNorm[tidx]);
           }


	}

	/* Add the integrand where appropriate */ 
	for ( int vol=0; vol<intvolumes; vol++)
	{
            CCTK_INT my_calc_norms_every = calc_norms_every_vol[vol] >= 0 ? calc_norms_every_vol[vol] : calc_norms_every;
            assert(my_calc_norms_every >= 0);
            if ( (my_calc_norms_every == 0) || (cctk_iteration % my_calc_norms_every != 0) )
                continue;

	    weightL = VolWeights[vol][idx];
            if ( weightL==0 || (weightL > 0 && omit_all_buffers[vol]>0) ) {
               continue;
            }

	    r2L = SQR(xL-x0L[vol]) + SQR(yL-y0L[vol]) + SQR(zL-z0L[vol]) ;
	    if ( r2L < SQR(rad0L[vol]) )
	    {
               for ( int nscalar=0; nscalar<num_scalars; nscalar++ ) {
                   localNorms[ vol*norms_per_vol + nscalar ] += dNorm[ nscalar ]*ExtraWeightVal[nscalar]*weightL*sqrt(detg)*dV; 
                   localCount[ vol*norms_per_vol + nscalar ] +=                  ExtraWeightVal[nscalar]*weightL*sqrt(detg)*dV;
               } 
               for ( int norm_in_vol=num_scalars; norm_in_vol<norms_per_vol; norm_in_vol++ ) {
                   localNorms[ vol*norms_per_vol + norm_in_vol ] += dNorm[ norm_in_vol ]*weightL*sqrt(detg)*dV; 
                   localCount[ vol*norms_per_vol + norm_in_vol ] +=                      weightL*sqrt(detg)*dV; 
               }
               if (verbose && (idx%100==0) ) {
                  CCTK_VInfo(CCTK_THORNSTRING,"Adding dN at r=%g (<%g) for vol %d. weight=%g, sqrt(detg)=%g, dV=%g",
                        sqrt(r2L),rad0L[vol],vol,weightL,sqrt(detg),dV);
               }
               if (debug && (idx%100==0) ) {
                  CCTK_VInfo(CCTK_THORNSTRING,"Adding dN at r=%g (<%g) for vol %d. weight=%g",sqrt(r2L),rad0L[vol],vol,weightL);
                  CCTK_VInfo(CCTK_THORNSTRING," Current Local Scalar Norm: (vol, ||alp||)={(0,%g),(1,%g),(2,%g),(3,%g)}", localNorms[ 0 ],
                      localNorms[ norms_per_vol ], localNorms[ 2*norms_per_vol ], localNorms[ 3*norms_per_vol ]);
                  CCTK_VInfo(CCTK_THORNSTRING," Current Local Tensor Norm: (vol, ||g||)={(0,%g),(1,%g),(2,%g),(3,%g)}", localNorms[ num_scalars+num_vectors ],
                      localNorms[ norms_per_vol + (num_scalars+num_vectors) ], localNorms[ 2*norms_per_vol + (num_scalars+num_vectors) ], 
                      localNorms[ norms_per_vol + (num_scalars+num_vectors) ]);
               }
	    }
        }
      }
    }
  }

 if (debug) {
    CCTK_VInfo(CCTK_THORNSTRING,"Local Scalar Norm: (vol, ||alp||)={(0,%g),(1,%g),(2,%g),(3,%g)}", localNorms[ 0 ],
               localNorms[ norms_per_vol ], localNorms[ 2*norms_per_vol ], localNorms[ 3*norms_per_vol ]);
    CCTK_VInfo(CCTK_THORNSTRING,"Local Tensor Norm: (vol, ||g||)={(0,%g),(1,%g),(2,%g),(3,%g)}", localNorms[ num_scalars+num_vectors ],
               localNorms[ norms_per_vol + (num_scalars+num_vectors) ], localNorms[ 2*norms_per_vol + (num_scalars+num_vectors) ], 
               localNorms[ norms_per_vol + (num_scalars+num_vectors) ]);
 }
 return 0;

}

/* IO */
void SubspaceNorms_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char fmode[2];
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  char addformat[50];
  FILE *file;


  int did_calc_norm=0;
  int my_calc_norm_every;
  int calculated_new_norm[intvolumes];
  for ( int vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_norm_every = calc_norms_every_vol[vol] >= 0 ? calc_norms_every_vol[vol] : calc_norms_every;
    assert(my_calc_norm_every >= 0);
    if ( (my_calc_norm_every > 0) && (cctk_iteration % my_calc_norm_every == 0) ) {
        did_calc_norm = 1;
        calculated_new_norm[vol]=1;
    } else {
        calculated_new_norm[vol]=0;
    }
  }

  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( (myproc == 0) && did_calc_norm )
  {

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "Writing output");
    }

    CCTK_INT norms_per_vol = num_scalars + num_vectors + num_tensors;
    for ( int vol=0 ; vol<intvolumes; vol++ ) { 

       if ( calculated_new_norm[vol] == 0 ) {
          // Skipping this volume
          continue;
       }

       // One file per integration volume
       sprintf(varname, "SubspaceNorms_vol_%d", vol);

       filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                                strlen (file_extension) +2);
       assert(filename);
       sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

       // Open file, appending if already written.
       if ( file_created[vol]>0 ) {
          sprintf(fmode,"%s","a");
       } else {
          sprintf(fmode,"%s","w");
       }
       file = fopen (filename, fmode);
       if (!file) {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open norm output file '%s'",
                  filename);
       }

       // Write header on startup
       if (file_created[vol]<=0) {
         fprintf(file,"# Subspace Integrals, intvolume=%d\n",vol);
         fprintf(file,"# gnuplot column index:\n");
         fprintf(file,"# 1:it 2:t \n");

         fprintf(file,"#    Scalar Integrals: ");
         for ( int nscalar=0; nscalar<num_scalars; nscalar++) 
         {
             if ( norm_exponent[nscalar] == 2 ) {
	        fprintf(file," %d:%s", 3+nscalar, CCTK_VarName(varinfo_scalars[nscalar]));
             } else {
	        fprintf(file," %d:%s(N%d)", 3+nscalar, CCTK_VarName(varinfo_scalars[nscalar]), norm_exponent[nscalar]);
             }
         }
         fprintf(file,"\n");

         fprintf(file,"#    Vector Integrals: ");
         for ( int nvector=0; nvector<num_vectors; nvector++) 
         {
	     fprintf(file," %d:%s", 3+num_scalars+nvector, CCTK_GroupName(groupinfo_vectors[nvector]));
         }
         fprintf(file,"\n");

         fprintf(file,"#    Tensor Integrals: ");
         for ( int ntensor=0; ntensor<num_tensors; ntensor++) 
         {
	     fprintf(file," %d:%s", 3+num_scalars+num_vectors+ntensor, CCTK_GroupName(groupinfo_tensors[ntensor]));
         }
         fprintf(file,"\n");

         if ( output_ncount ) {

	    int coloffset=3+(num_scalars+num_vectors+num_tensors);
            fprintf(file,"#    Counters ( ie. Sum(w detg dV) ): ");
            fprintf(file,"##     Scalars: ");
            for ( int nscalar=0; nscalar<num_scalars; nscalar++) 
            {
	        fprintf(file," %d:V(%s)", coloffset+nscalar, CCTK_VarName(varinfo_scalars[nscalar]));
            }
            fprintf(file,"\n");
            fprintf(file,"##     Vectors: ");
            for ( int nvector=0; nvector<num_vectors; nvector++) 
            {
  	        fprintf(file," %d:V(%s)", coloffset+num_scalars+nvector, CCTK_GroupName(groupinfo_vectors[nvector]));
            }
            fprintf(file,"\n");
            fprintf(file,"##     Tensors: ");
            for ( int ntensor=0; ntensor<num_tensors; ntensor++) 
            {
	        fprintf(file," %d:V(%s)", coloffset+num_scalars+num_vectors+ntensor, CCTK_GroupName(groupinfo_tensors[ntensor]));
            }
            fprintf(file,"\n");

         }

       } // End header information

       // Format string for data output
       sprintf (format_str_real, "%%d\t%%%s",out_format);
       sprintf (addformat,"\t%%%s",out_format);

       // Write data to file
       fprintf(file, format_str_real, cctk_iteration, cctk_time);
       for ( int nscalar=0; nscalar<num_scalars; nscalar++) {
           CCTK_INT nidx = vol*norms_per_vol + nscalar;
           CCTK_REAL outdata = ( norm_as_output > 0 ? ScalarNorms[vol*num_scalars+nscalar]/NormCounter[nidx] : ScalarNorms[vol*num_scalars+nscalar] );
           if ( norm_exponent[nscalar] == 2 ) {
              fprintf(file,addformat,sqrt(outdata));
           } else {
              fprintf(file,addformat,outdata);
           }
       }
       for ( int nvector=0; nvector<num_vectors; nvector++) {
           CCTK_INT nidx = vol*norms_per_vol + num_scalars + nvector;
           CCTK_REAL outdata = ( norm_as_output > 0 ? VectorNorms[vol*num_vectors+nvector]/NormCounter[nidx] : VectorNorms[vol*num_vectors+nvector] );
           fprintf(file,addformat, sqrt(outdata));
       }
       for ( int ntensor=0; ntensor<num_tensors; ntensor++) {
           CCTK_INT nidx = vol*norms_per_vol + num_scalars + num_vectors + ntensor;
           CCTK_REAL outdata = ( norm_as_output > 0 ? TensorNorms[vol*num_tensors+ntensor]/NormCounter[nidx] : TensorNorms[vol*num_tensors+ntensor] );
           fprintf(file,addformat, sqrt(outdata));
       }
       if ( output_ncount ) {
           for ( int nscalar=0; nscalar<num_scalars; nscalar++) {
               fprintf(file,addformat,NormCounter[vol*norms_per_vol + nscalar]);
           }
           for ( int nvector=0; nvector<num_vectors; nvector++) {
               fprintf(file,addformat,NormCounter[vol*norms_per_vol + num_scalars + nvector]);
           }
           for ( int ntensor=0; ntensor<num_tensors; ntensor++) {
               fprintf(file,addformat,NormCounter[vol*norms_per_vol + num_scalars + num_vectors + ntensor]);
           }
       }
       fprintf(file,"\n");

       fclose(file);
       free(filename);

       if (file_created[vol]==0) {
          file_created[vol]=1;
       }

    }
  }
}
