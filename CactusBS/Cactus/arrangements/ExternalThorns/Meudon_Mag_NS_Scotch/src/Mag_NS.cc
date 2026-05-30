#include <cassert>
#include <cstdio>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mag_ns.h>

#define TINY 1e-15
#define SQR(x) ((x)*(x))

using namespace std;



extern "C"
void ID_Mag_NS_Scotch_initialise (CCTK_ARGUMENTS)
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
  CCTK_REAL const mB      = 1.66e-27; // unit baryon mass [kg]
  CCTK_REAL const rho_nuc = mB*1.e44; // nuclear density mB*n_nuc [kg/m^3]
  CCTK_REAL const j_Lorene_unit  = 1e11; // Units of A/m^2
  CCTK_REAL const Avec_fudge_trans = 2*M_PI*1.078; // FIX: Don't know why.

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
  //CCTK_REAL const A_unit     = B_unit*cactusL/( mag_unit*r_Lorene_unit );
  CCTK_REAL const A_unit     = B_unit*cactusL/(Avec_fudge_trans); //( r_Lorene_unit );
  CCTK_REAL       K_unit     = 0; // Needs gamma from star to calculate
  
  // CCTK_EQUALS is expensive
  CCTK_INT const set_Bvec = CCTK_EQUALS( initial_Bvec, "ID_Mag_NS_Scotch" ); 
  CCTK_INT const set_Avec = CCTK_EQUALS( initial_Avec, "ID_Mag_NS_Scotch" ); 
  CCTK_INT const set_Aphi = CCTK_EQUALS( initial_Aphi, "ID_Mag_NS_Scotch" ); 

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

  if ( mag_ns.gamma_poly > 0 ) {
       // Lorene units of [rho_nuc^(1-gamma) c^2]_SI
       K_unit = mag_ns.kappa_poly * pow( (rho_nuc/rho_unit) , 1.-mag_ns.gamma_poly );
       CCTK_VInfo( CCTK_THORNSTRING, "Polytrope, kappa_poly = %g in Cactus units", K_unit );
  }
  
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
    
    eps[i] = mag_ns.ener_spec[i];
    //eps[i] = rho[i] * mag_ns.ener_spec[i] / ener_unit;
    
    velx[i] = mag_ns.u_euler_x[i] / vel_unit;
    vely[i] = mag_ns.u_euler_y[i] / vel_unit;
    velz[i] = mag_ns.u_euler_z[i] / vel_unit;
    
    if ( set_Bvec > 0 ) {
       Bvec[i          ] = mag_ns.bb_x[i] / B_unit;
       Bvec[i+  npoints] = mag_ns.bb_y[i] / B_unit;
       Bvec[i+2*npoints] = mag_ns.bb_z[i] / B_unit;
    }

    if ( set_Avec > 0 ) {
       double Rcyl2 = x[i]*x[i] + y[i]*y[i] + TINY;
       Avec[i          ] = -(y[i]*mag_ns.Aphi[i])/(Rcyl2*A_unit);
       Avec[i+  npoints] = x[i]*mag_ns.Aphi[i]/(Rcyl2*A_unit);
       Avec[i+2*npoints] = 0;;
       //Avec[i          ] = mag_ns.Av_x[i]/A_unit; 
       //Avec[i+  npoints] = mag_ns.Av_y[i]/A_unit;
       //Avec[i+2*npoints] = mag_ns.Av_z[i]/A_unit;
    }

    if ( set_Aphi > 0 ) {
       Aphi[i] = mag_ns.At[i]/A_unit;
    }

    /****************************/
    /* Fill auxiliary variables */

    // w_lorentz filled in during initial prim 2 con
    // Initial divclean_psi as well?

    if ( divclean_psi != NULL ) {
       divclean_psi[i] = 0.;
    }

    if ( mag_ns.gamma_poly > 0 ) {

       press[i] = K_unit * pow( rho[i], mag_ns.gamma_poly );
       //press[i] = rho[i] * eps[i] * ( mag_ns.gamma_poly - 1. );

       /*if ( rho[i] > 0.001 ) {
          CCTK_VInfo(CCTK_THORNSTRING, "COMPARE: rho=%g, eps=%g, P = %g (rho*eps*(G-1)) or %g (K*rho^G) with K=%g",
                     rho[i], eps[i], press[i],  rho[i] * eps[i] * ( mag_ns.gamma_poly - 1. ), K_unit);
       }*/

    } else {
       CCTK_WARN(CCTK_WARN_ABORT, "Not a polytrope LORENE star! Add pressure calculation!"); 
    }
    /****************************/

    /* Impose atmosphere */
    if ( rho[i] < *whisky_rho_min ) {
       if ( *whisky_rho_min > 10. ) {
          CCTK_WARN(CCTK_WARN_ABORT, "Atmosphere isn't properly set yet!");
       }
       rho[i] = *whisky_rho_min;
       velx[i] = 0.;
       vely[i] = 0.;
       velz[i] = 0.;
       press[i] = K_unit * pow( *whisky_rho_min, mag_ns.gamma_poly );
       eps[i] = press[i]/( *whisky_rho_min * ( mag_ns.gamma_poly - 1.) );
    }


  } // for i
  
  
  
  CCTK_INFO ("Setting time derivatives of lapse and shift");
  {
    // These initial data assume stationarity
    
    if (CCTK_EQUALS (initial_dtlapse, "ID_Mag_NS_Scotch")) {
#pragma omp parallel for
      for (int i=0; i<npoints; ++i) {
        dtalp[i] = 0.0;
      }
    } else if (CCTK_EQUALS (initial_dtlapse, "none")) {
      // do nothing
    } else {
      CCTK_WARN (CCTK_WARN_ABORT, "internal error");
    }
    
    if (CCTK_EQUALS (initial_dtshift, "ID_Mag_NS_Scotch")) {
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

