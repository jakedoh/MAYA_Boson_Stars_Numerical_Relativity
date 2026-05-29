#include <assert.h>

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 
#include "cctk_Functions.h" 

#define DEBUG

#include "SpaceMask.h"

namespace AHMask {

const int max_centres = 3;


void get_centre_parameters(CCTK_ARGUMENTS, bool initial, int c, 
                                  CCTK_REAL *radius, CCTK_REAL coords[3]);

// some math helpers
static inline double sqr(const double x) {return x*x;}

// the main callable function here, set weigths to the intersection of
// CarpetReduce::weights and is_outside_of_the_centres
static void Internal_Set_Weight(CCTK_ARGUMENTS, bool first_call)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;    

  CCTK_REAL radii[max_centres], centre_coords[max_centres][3];
  CCTK_REAL xx,yy,zz,rr[max_centres];
  CCTK_REAL *CarpetWeights;
  CCTK_INT ah_bits=0, outside=0, inside=0, *spacemask=NULL; // pacify gcc
  bool is_not_outside, is_inside;

  // If we use Carpet then obtain a pointer to rhe carpet variable. This will be CarpetWeights
  CarpetWeights = NULL;
  if(CCTK_IsThornActive("CarpetReduce")!=0)
  {
    CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,"CarpetReduce::weight"));
    assert(CarpetWeights);
  }

  // are we supposed to incorporate the ExcisionMask into this?
  if(use_space_mask)
  {
    ah_bits = SpaceMask_GetTypeBits(mask_bitfield_name);
    inside = SpaceMask_GetStateBits(mask_bitfield_name, space_mask_inside_value);
    outside = SpaceMask_GetStateBits(mask_bitfield_name, space_mask_outside_value);
    spacemask = (CCTK_INT *)CCTK_VarDataPtr(cctkGH,0,space_mask_gridfn_name);
    if(spacemask == NULL)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Cannot find grid "
              "function '%s'. Did you turn on storage for it?",
              space_mask_gridfn_name);
      /* NOTREACHED */
    }
    assert(spacemask != NULL);
    assert(outside != -1);
    assert(ah_bits != 0);
  }

  for(int c = 0 ; c < number_of_centres ; c++)
  {
    // coordinates of centre in question
    get_centre_parameters(CCTK_PASS_CTOC, first_call, c, &radii[c], centre_coords[c]);

    // save this for later speedup
    rr[c] = sqr(radii[c]);
  }

  long int inside_points = 0;
  long int not_outside_points = 0;
  for(int k = 0 ; k < cctk_lsh[2] ; k++)
  {
    for(int j = 0 ; j < cctk_lsh[1] ; j++)
    {
      for(int i = 0 ; i < cctk_lsh[0] ; i++)
      {
        const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);

        is_not_outside = false;
        is_inside = false;

        // check if we are inside any of the centres
        for(int c = 0  ; c < number_of_centres ; c++)
        {
          xx = sqr(x[idx] - centre_coords[c][0]);
          yy = sqr(y[idx] - centre_coords[c][1]);
          zz = sqr(z[idx] - centre_coords[c][2]);

          if(xx+yy+zz < rr[c]) {
            is_not_outside = true;
            is_inside = true;
            break;
          }
        }

        // check if we are inside of the excised region
        if(use_space_mask && !SpaceMask_CheckStateBits(spacemask, idx, ah_bits, outside))
        {
            ++not_outside_points;
            is_not_outside = true;
        }
        if(use_space_mask && SpaceMask_CheckStateBits(spacemask, idx, ah_bits, inside))
        {
            ++inside_points;
            is_inside = true;
        }

        // if we are inside take note of this, else use Carpet's weight
        if(is_inside)
          AHmasked_weight[idx] = 0.;
        else if(CarpetWeights != NULL)
          AHmasked_weight[idx] = CarpetWeights[idx];
        else 
          AHmasked_weight[idx] = 1.;

        if(create_buffered_mask)
        {
          if(is_not_outside)
            AHmask_buffered[idx] = 0.;
          else if(CarpetWeights != NULL)
            AHmask_buffered[idx] = CarpetWeights[idx];
          else
            AHmask_buffered[idx] = 1.;
        }

      }
    }
  }
  if(verbose > 1)
    CCTK_VInfo(CCTK_THORNSTRING, "inside_points: %ld, buffered_points: %ld", inside_points, not_outside_points);
}

// until the spherical surface are set for the first time we have to use the
// coodinates given in the the paramters (same as RegridBoxes)
extern "C" void AHMask_Set_Weight_Initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  Internal_Set_Weight(CCTK_PASS_CTOC, true);
}

extern "C" void AHMask_Set_Weight(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  Internal_Set_Weight(CCTK_PASS_CTOC, false);
}
}
