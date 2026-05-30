/* (c) 2009 Erik Schnetter
 * (c) 2010 Frank Loeffler
 *
 * Adapted for Scotch by Christopher W. Evans */

#include <cassert>
#include <cstdio>
#include <ios>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "EOS_Base.h"
#include <bin_ns.h>
#include <unites.h>

using namespace std;

static void set_dt_from_domega(CCTK_ARGUMENTS, CCTK_REAL const *const var,
                               CCTK_REAL *const dtvar, CCTK_REAL const &omega) {
  DECLARE_CCTK_ARGUMENTS;

  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<CCTK_REAL> dxvar(npoints), dyvar(npoints);

  Diff_gv(cctkGH, 0, var, &dxvar[0], -1);
  Diff_gv(cctkGH, 1, var, &dyvar[0], -1);

#pragma omp parallel for
  for (int i = 0; i < npoints; ++i) {
    CCTK_REAL const ephix = +y[i];
    CCTK_REAL const ephiy = -x[i];
    CCTK_REAL const dphi_var = ephix * dxvar[i] + ephiy * dyvar[i];
    dtvar[i] = omega * dphi_var;
  }
}

extern "C" void Meudon_Bin_NS_Scotch_initialise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Setting up LORENE Bin_NS initial data");

  // Meudon data are distributed in SI units (MKSA).  Here are some
  // conversion factors.
  // Be aware: these are the constants Lorene uses. They do differ from other
  // conventions, but they gave the best results in some tests.

  CCTK_REAL const c_light = Unites::c_si; // speed of light [m/s]
  CCTK_REAL const nuc_dens =
      Unites::rhonuc_si; // Nuclear density as used in Lorene units [kg/m^3]
  CCTK_REAL const G_grav = Unites::g_si; // gravitational constant [m^3/kg/s^2]
  CCTK_REAL const M_sun = Unites::msol_si; // solar mass [kg]

  // Cactus units in terms of SI units:
  // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
  // for the magnetic field)
  CCTK_REAL const cactusM = M_sun * units;
  CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light, 2);
  CCTK_REAL const cactusT = cactusL / c_light;

  // Other quantities in terms of Cactus units
  CCTK_REAL const coord_unit = cactusL / 1.0e+3;        // from km (~1.477)
  CCTK_REAL const rho_unit = cactusM / pow(cactusL, 3); // from kg/m^3
  CCTK_REAL const rho_atmosphere_code = rho_atmosphere * units * units;

  CCTK_INFO("Setting up coordinates");

  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<double> xx(npoints), yy(npoints), zz(npoints);

#pragma omp parallel for
  for (int i = 0; i < npoints; ++i) {
    xx[i] = x[i] * coord_unit;
    yy[i] = y[i] * coord_unit;
    zz[i] = z[i] * coord_unit;
  }

  CCTK_VInfo(CCTK_THORNSTRING, "Reading from file \"%s\"", filename);

  try {
    Bin_NS bin_ns(npoints, &xx[0], &yy[0], &zz[0], filename);

    CCTK_VInfo(CCTK_THORNSTRING, "omega [rad/s]:       %g", bin_ns.omega);
    CCTK_VInfo(CCTK_THORNSTRING, "dist [km]:           %g", bin_ns.dist);
    CCTK_VInfo(CCTK_THORNSTRING, "dist_mass [km]:      %g", bin_ns.dist_mass);
    CCTK_VInfo(CCTK_THORNSTRING, "mass1_b [M_sun]:     %g", bin_ns.mass1_b);
    CCTK_VInfo(CCTK_THORNSTRING, "mass2_b [M_sun]:     %g", bin_ns.mass2_b);
    CCTK_VInfo(CCTK_THORNSTRING, "mass_ADM [M_sun]:    %g", bin_ns.mass_adm);
    CCTK_VInfo(CCTK_THORNSTRING, "L_tot [G M_sun^2/c]: %g", bin_ns.angu_mom);
    CCTK_VInfo(CCTK_THORNSTRING, "rad1_x_comp [km]:    %g", bin_ns.rad1_x_comp);
    CCTK_VInfo(CCTK_THORNSTRING, "rad1_y [km]:         %g", bin_ns.rad1_y);
    CCTK_VInfo(CCTK_THORNSTRING, "rad1_z [km]:         %g", bin_ns.rad1_z);
    CCTK_VInfo(CCTK_THORNSTRING, "rad1_x_opp [km]:     %g", bin_ns.rad1_x_opp);
    CCTK_VInfo(CCTK_THORNSTRING, "rad2_x_comp [km]:    %g", bin_ns.rad2_x_comp);
    CCTK_VInfo(CCTK_THORNSTRING, "rad2_y [km]:         %g", bin_ns.rad2_y);
    CCTK_VInfo(CCTK_THORNSTRING, "rad2_z [km]:         %g", bin_ns.rad2_z);
    CCTK_VInfo(CCTK_THORNSTRING, "rad2_x_opp [km]:     %g", bin_ns.rad2_x_opp);

    assert(bin_ns.np == npoints);

    CCTK_INFO("Filling in Cactus grid points");

#pragma omp parallel for
    for (int i = 0; i < npoints; ++i) {

      if (CCTK_EQUALS(initial_lapse, "Meudon_Bin_NS_Scotch")) {
        alp[i] = bin_ns.nnn[i];
      }

      if (CCTK_EQUALS(initial_shift, "Meudon_Bin_NS_Scotch")) {
        betax[i] = -bin_ns.beta_x[i];
        betay[i] = -bin_ns.beta_y[i];
        betaz[i] = -bin_ns.beta_z[i];
      }

      if (CCTK_EQUALS(initial_data, "Meudon_Bin_NS_Scotch")) {
        gxx[i] = bin_ns.g_xx[i];
        gxy[i] = bin_ns.g_xy[i];
        gxz[i] = bin_ns.g_xz[i];
        gyy[i] = bin_ns.g_yy[i];
        gyz[i] = bin_ns.g_yz[i];
        gzz[i] = bin_ns.g_zz[i];

        kxx[i] = bin_ns.k_xx[i] * coord_unit;
        kxy[i] = bin_ns.k_xy[i] * coord_unit;
        kxz[i] = bin_ns.k_xz[i] * coord_unit;
        kyy[i] = bin_ns.k_yy[i] * coord_unit;
        kyz[i] = bin_ns.k_yz[i] * coord_unit;
        kzz[i] = bin_ns.k_zz[i] * coord_unit;
      }

      if (CCTK_EQUALS(initial_data, "Meudon_Bin_NS_Scotch")) {
        rho[i] = bin_ns.nbar[i] / rho_unit;
        eps[i] = bin_ns.ener_spec[i];
        press[i] = EOS_Pressure(*whisky_eos_handle, rho[i], eps[i]);

        velx[i] = bin_ns.u_euler_x[i];
        vely[i] = bin_ns.u_euler_y[i];
        velz[i] = bin_ns.u_euler_z[i];

        // Especially the velocity is set to strange values outside of the
        // matter region, so take care of this in the following way
        if (rho[i] < rho_atmosphere_code) {
          rho[i] = rho_atmosphere_code;
          velx[i] = 0.0;
          vely[i] = 0.0;
          velz[i] = 0.0;
          press[i] =
              EOS_Pressure(*whisky_polytrope_handle, rho_atmosphere_code, 0.);
          eps[i] = EOS_SpecificIntEnergy(*whisky_polytrope_handle,
                                         rho_atmosphere_code, press[i]);
        }
      }
    }

    {
      // Angular velocity
      CCTK_REAL const omega = bin_ns.omega * cactusT;

      // These initial data assume a helical Killing vector field

      if (CCTK_EQUALS(initial_lapse, "Meudon_Bin_NS_Scotch")) {
        if (CCTK_EQUALS(initial_dtlapse, "Meudon_Bin_NS_Scotch")) {
          CCTK_INFO("Calculating time derivatives of lapse");
          set_dt_from_domega(CCTK_PASS_CTOC, alp, dtalp, omega);
        } else if (CCTK_EQUALS(initial_dtlapse, "none")) {
          // do nothing
        } else {
          CCTK_WARN(CCTK_WARN_ABORT, "internal error");
        }
      }

      if (CCTK_EQUALS(initial_shift, "Meudon_Bin_NS_Scotch")) {
        if (CCTK_EQUALS(initial_dtshift, "Meudon_Bin_NS_Scotch")) {
          CCTK_INFO("Calculating time derivatives of shift");
          set_dt_from_domega(CCTK_PASS_CTOC, betax, dtbetax, omega);
          set_dt_from_domega(CCTK_PASS_CTOC, betay, dtbetay, omega);
          set_dt_from_domega(CCTK_PASS_CTOC, betaz, dtbetaz, omega);
        } else if (CCTK_EQUALS(initial_dtshift, "none")) {
          // do nothing
        } else {
          CCTK_WARN(CCTK_WARN_ABORT, "internal error");
        }
      }
    }

    CCTK_INFO("Done.");
  } catch (ios::failure e) {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not read initial data from file '%s': %s", filename,
               e.what());
  }
}
