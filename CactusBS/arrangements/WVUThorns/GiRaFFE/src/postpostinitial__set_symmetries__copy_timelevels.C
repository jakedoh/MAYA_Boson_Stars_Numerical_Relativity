//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "GiRaFFE_headers.h"

extern "C" void GiRaFFE_PostPostInitial_Set_Symmetries__Copy_Timelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  //For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  // or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE AND PRIMIIVE VARIABLES!
    int ierr;
    ierr=CartSymGN(cctkGH,"GiRaFFE::grmhd_conservatives"); if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"Microsoft error code #1874109358120048. Grep it in the source code");
    ierr=CartSymGN(cctkGH,"GiRaFFE::grmhd_primitives_allbutBi"); if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"Microsoft error code #1874109358120049. Grep it in the source code");

    // Finish up by setting symmetry ghostzones on Bx, By, Bz, and their staggered variants.
    const CCTK_REAL gridfunc_syms_Bx[3] = {-1, 1,-Sym_Bz};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bx        , gridfunc_syms_Bx,0,0,0);
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bx_stagger, gridfunc_syms_Bx,1,0,0);
    const CCTK_REAL gridfunc_syms_By[3] = { 1,-1,-Sym_Bz};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  By        , gridfunc_syms_Bx,0,0,0);
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  By_stagger, gridfunc_syms_By,0,1,0);
    const CCTK_REAL gridfunc_syms_Bz[3] = { 1, 1, Sym_Bz};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bz        , gridfunc_syms_Bz,0,0,0);
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bz_stagger, gridfunc_syms_Bz,0,0,1);

    const CCTK_REAL gridfunc_syms_psi6phi[3] = { 1, 1,      1};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,psi6phi     , gridfunc_syms_psi6phi,1,1,1);
    const CCTK_REAL gridfunc_syms_Ax[3]      = {-1, 1, Sym_Bz};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Ax        , gridfunc_syms_Ax,0,1,1);
    const CCTK_REAL gridfunc_syms_Ay[3]      = { 1,-1, Sym_Bz};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Ay        , gridfunc_syms_Ay,1,0,1);
    const CCTK_REAL gridfunc_syms_Az[3]      = { 1, 1,-Sym_Bz};
    GiRaFFE_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Az        , gridfunc_syms_Az,1,1,0);
  }


  //------------------------------------------------------------------
  // FILL _p AND _p_p TIMELEVELS. Probably don't need to do this if
  // Carpet::init_fill_timelevels=yes  and
  // MoL::initial_data_is_crap = yes
  // NOTE: We don't fill metric data here.
  // FIXME: Do we really need this?

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        mhd_st_x_p[index] = mhd_st_x[index];
        mhd_st_y_p[index] = mhd_st_y[index];
        mhd_st_z_p[index] = mhd_st_z[index];

        psi6phi_p[index] = psi6phi[index];
        Ax_p[index] = Ax[index];
        Ay_p[index] = Ay[index];
        Az_p[index] = Az[index];

        mhd_st_x_p_p[index] = mhd_st_x[index];
        mhd_st_y_p_p[index] = mhd_st_y[index];
        mhd_st_z_p_p[index] = mhd_st_z[index];

        psi6phi_p_p[index] = psi6phi[index];
        Ax_p_p[index] = Ax[index];
        Ay_p_p[index] = Ay[index];
        Az_p_p[index] = Az[index];
      }
}
