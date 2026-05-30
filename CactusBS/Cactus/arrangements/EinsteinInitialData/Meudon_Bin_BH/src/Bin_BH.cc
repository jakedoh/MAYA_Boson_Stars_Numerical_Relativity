#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <bin_bh.h>

using namespace std;


// define namespace here for old versions of Lorene that don't do so
namespace Lorene {}
using namespace Lorene;

static void set_dt_from_domega (CCTK_ARGUMENTS,
                                CCTK_REAL const* const var,
                                CCTK_REAL      * const dtvar,
                                CCTK_REAL const& omega)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<CCTK_REAL> dxvar(npoints), dyvar(npoints);
  
  Diff_gv (cctkGH, 0, var, &dxvar[0], -1);
  Diff_gv (cctkGH, 1, var, &dyvar[0], -1);
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    CCTK_REAL const ephix = +y[i];
    CCTK_REAL const ephiy = -x[i];
    CCTK_REAL const dphi_var = ephix * dxvar[i] + ephiy * dyvar[i];
    dtvar[i] = omega * dphi_var;
  }
}



extern "C"
void ID_Bin_BH_initialise (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INFO ("Setting up LORENE Bin_BH initial data");
  
  
  
  CCTK_INFO ("Setting up coordinates");
  
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<double> xx(npoints), yy(npoints), zz(npoints);
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    xx[i] = x[i];
    yy[i] = y[i];
    zz[i] = z[i];
  }
  
  
  
  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);
  
  try {
  Bin_BH bin_bh (npoints, &xx[0], &yy[0], &zz[0], 1, filename);
  
  CCTK_VInfo (CCTK_THORNSTRING, "Omega [1/a]: %g", bin_bh.omega);
  CCTK_VInfo (CCTK_THORNSTRING, "dist [a]:    %g", bin_bh.dist);
  CCTK_VInfo (CCTK_THORNSTRING, "radius1 [a]: %g", 1.0);
  CCTK_VInfo (CCTK_THORNSTRING, "radius2 [a]: %g", bin_bh.radius2);
  assert (bin_bh.np == npoints);
  
  
  
  CCTK_INFO ("Filling in Cactus grid points");
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    
    alp[i] = bin_bh.nnn[i];
    
    betax[i] = bin_bh.beta_x[i];
    betay[i] = bin_bh.beta_y[i];
    betaz[i] = bin_bh.beta_z[i];
    
    CCTK_REAL g[3][3];
    g[0][0] = bin_bh.g_xx[i];
    g[0][1] = bin_bh.g_xy[i];
    g[0][2] = bin_bh.g_xz[i];
    g[1][1] = bin_bh.g_yy[i];
    g[1][2] = bin_bh.g_yz[i];
    g[2][2] = bin_bh.g_zz[i];
    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    
    CCTK_REAL ku[3][3];
    ku[0][0] = bin_bh.k_xx[i];
    ku[0][1] = bin_bh.k_xy[i];
    ku[0][2] = bin_bh.k_xz[i];
    ku[1][1] = bin_bh.k_yy[i];
    ku[1][2] = bin_bh.k_yz[i];
    ku[2][2] = bin_bh.k_zz[i];
    ku[1][0] = ku[0][1];
    ku[2][0] = ku[0][2];
    ku[2][1] = ku[1][2];
    
    CCTK_REAL k[3][3];
    for (int a=0; a<3; ++a) {
      for (int b=0; b<3; ++b) {
        k[a][b] = 0.0;
        for (int c=0; c<3; ++c) {
          for (int d=0; d<3; ++d) {
            k[a][b] += g[a][c] * g[b][d] * ku[c][d];
          }
        }
      }
    }
    
    gxx[i] = g[0][0];
    gxy[i] = g[0][1];
    gxz[i] = g[0][2];
    gyy[i] = g[1][1];
    gyz[i] = g[1][2];
    gzz[i] = g[2][2];
    
    kxx[i] = k[0][0];
    kxy[i] = k[0][1];
    kxz[i] = k[0][2];
    kyy[i] = k[1][1];
    kyz[i] = k[1][2];
    kzz[i] = k[2][2];
    
  } // for i
  
  
  
  CCTK_INFO ("Calculating time derivatives of lapse and shift");
  {
    // Angular velocity
    CCTK_REAL const omega = bin_bh.omega * bin_bh.dist;
    
    // These initial data assume a helical Killing vector field
    
    if (CCTK_EQUALS (initial_dtlapse, "ID_Bin_BH")) {
      set_dt_from_domega (CCTK_PASS_CTOC, alp, dtalp, omega);
    } else if (CCTK_EQUALS (initial_dtlapse, "none")) {
      // do nothing
    } else {
      CCTK_WARN (CCTK_WARN_ABORT, "internal error");
    }
    
    if (CCTK_EQUALS (initial_dtshift, "ID_Bin_BH")) {
      set_dt_from_domega (CCTK_PASS_CTOC, betax, dtbetax, omega);
      set_dt_from_domega (CCTK_PASS_CTOC, betay, dtbetay, omega);
      set_dt_from_domega (CCTK_PASS_CTOC, betaz, dtbetaz, omega);
    } else if (CCTK_EQUALS (initial_dtshift, "none")) {
      // do nothing
    } else {
      CCTK_WARN (CCTK_WARN_ABORT, "internal error");
    }
  }
  
  
  
  CCTK_INFO ("Done.");
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial data from file '%s': %s", filename, e.what());
  }
}
