#include <cassert>
#include <cstdio>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mag_ns.h>

using namespace std;


// define namespace here for old versions of Lorene that don't do so
namespace Lorene {}
using namespace Lorene;

extern "C"
void ID_Mag_NS_initialise (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INFO ("Setting up LORENE Mag_NS initial data");
  
  
  
  // Meudon data are distributed in SI units (MKSA).  Here are some
  // conversion factors.
  
  // Defined constants
  CCTK_REAL const c_light = 299792458.0; // speed of light [m/s]
  CCTK_REAL const mu0     = 4*M_PI * 1.0e-7; // vacuum permeability [N/A^2]
  CCTK_REAL const eps0    = 1 / (mu0 * pow(c_light,2));
  
  // Constants of nature (IAU, CODATA):
  CCTK_REAL const G_grav = 6.67428e-11; // gravitational constant [m^3/kg/s^2]
  CCTK_REAL const M_sun  = 1.98892e+30; // solar mass [kg]
  
  // Cactus units in terms of SI units:
  // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
  // for the magnetic field)
  CCTK_REAL const cactusM = M_sun;
  CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light,2);
  CCTK_REAL const cactusT = cactusL / c_light;
  CCTK_REAL const cactusB =
    sqrt(4*M_PI) / cactusL / sqrt(4*M_PI * eps0 * G_grav / pow(c_light,2));
  
  // Other quantities in terms of Cactus units
  CCTK_REAL const coord_unit = cactusL / 1.0e+3;         // from km
  CCTK_REAL const rho_unit   = cactusM / pow(cactusL,3); // from kg/m^3
  CCTK_REAL const ener_unit  = pow(cactusL,2);           // from c^2
  CCTK_REAL const vel_unit   = cactusL / cactusT / c_light; // from c
  CCTK_REAL const B_unit     = cactusB / 1.0e+9; // from 10^9 T
  
  
  
  CCTK_INFO ("Setting up coordinates");
  
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  
  vector<double> xx(npoints), yy(npoints), zz(npoints);
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    xx[i] = x[i] * coord_unit;
    yy[i] = y[i] * coord_unit;
    zz[i] = z[i] * coord_unit;
  }
  
  
  
  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);
  
  try {
  Mag_NS mag_ns (npoints, &xx[0], &yy[0], &zz[0], filename);
  
  CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:     %g", mag_ns.omega);
  CCTK_VInfo (CCTK_THORNSTRING, "rho_c [kg/m^3]     %g", mag_ns.rho_c);
  CCTK_VInfo (CCTK_THORNSTRING, "eps_c [c^2]:       %g", mag_ns.eps_c);
  CCTK_VInfo (CCTK_THORNSTRING, "mass_b [M_sun]:    %g", mag_ns.mass_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass_g [M_sun]:    %g", mag_ns.mass_g);
  CCTK_VInfo (CCTK_THORNSTRING, "r_eq [km]:         %g", mag_ns.r_eq);
  CCTK_VInfo (CCTK_THORNSTRING, "r_p [km]:          %g", mag_ns.r_p);
  CCTK_VInfo (CCTK_THORNSTRING, "L [G M_sum^2/c]:   %g", mag_ns.angu_mom);
  CCTK_VInfo (CCTK_THORNSTRING, "T/W:               %g", mag_ns.T_over_W);
  CCTK_VInfo (CCTK_THORNSTRING, "mu [A m^2]:        %g", mag_ns.magn_mom);
  CCTK_VInfo (CCTK_THORNSTRING, "B_z_pole [10^9 T]: %g", mag_ns.b_z_pole);
  CCTK_VInfo (CCTK_THORNSTRING, "B_z_eq [10^9 T]:   %g", mag_ns.b_z_eq);
  assert (mag_ns.np == npoints);
  
  
  
  CCTK_INFO ("Filling in Cactus grid points");
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    
    alp[i] = mag_ns.nnn[i];
    
    betax[i] = mag_ns.beta_x[i];
    betay[i] = mag_ns.beta_y[i];
    betaz[i] = mag_ns.beta_z[i];
    
    CCTK_REAL g[3][3];
    g[0][0] = mag_ns.g_xx[i];
    g[0][1] = mag_ns.g_xy[i];
    g[0][2] = mag_ns.g_xz[i];
    g[1][1] = mag_ns.g_yy[i];
    g[1][2] = mag_ns.g_yz[i];
    g[2][2] = mag_ns.g_zz[i];
    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    
    CCTK_REAL ku[3][3];
    ku[0][0] = mag_ns.k_xx[i];
    ku[0][1] = mag_ns.k_xy[i];
    ku[0][2] = mag_ns.k_xz[i];
    ku[1][1] = mag_ns.k_yy[i];
    ku[1][2] = mag_ns.k_yz[i];
    ku[2][2] = mag_ns.k_zz[i];
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
    
    rho[i] = mag_ns.nbar[i] / rho_unit;
    
    eps[i] = rho[i] * mag_ns.ener_spec[i] / ener_unit;
    
    vel[i          ] = mag_ns.u_euler_x[i] / vel_unit;
    vel[i+  npoints] = mag_ns.u_euler_y[i] / vel_unit;
    vel[i+2*npoints] = mag_ns.u_euler_z[i] / vel_unit;
    
    Bvec[i          ] = mag_ns.bb_x[i] / B_unit;
    Bvec[i+  npoints] = mag_ns.bb_y[i] / B_unit;
    Bvec[i+2*npoints] = mag_ns.bb_z[i] / B_unit;
    
  } // for i
  
  
  
  CCTK_INFO ("Setting time derivatives of lapse and shift");
  {
    // These initial data assume stationarity
    
    if (CCTK_EQUALS (initial_dtlapse, "ID_Mag_NS")) {
#pragma omp parallel for
      for (int i=0; i<npoints; ++i) {
        dtalp[i] = 0.0;
      }
    } else if (CCTK_EQUALS (initial_dtlapse, "none")) {
      // do nothing
    } else {
      CCTK_WARN (CCTK_WARN_ABORT, "internal error");
    }
    
    if (CCTK_EQUALS (initial_dtshift, "ID_Mag_NS")) {
#pragma omp parallel for
      for (int i=0; i<npoints; ++i) {
        dtbetax[i] = 0.0;
        dtbetay[i] = 0.0;
        dtbetaz[i] = 0.0;
      }
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
