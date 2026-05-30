
/*
  -------------------------------------------------------------------------------
  Copyright 2005 Scott C. Noble, Charles F. Gammie,
  Jonathan C. McKinney, and Luca Del Zanna


  This file is part of PVS-GRMHD.

  PVS-GRMHD is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PVS-GRMHD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PVS-GRMHD; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------------------
*/

// Function prototypes for this file:
static void raise_g(CCTK_REAL vcov[], CCTK_REAL gcon[][NDIM], CCTK_REAL vcon[]);
static void lower_g(CCTK_REAL vcon[], CCTK_REAL gcov[][NDIM], CCTK_REAL vcov[]);
static void ncov_calc(CCTK_REAL gcon[][NDIM],CCTK_REAL ncov[]) ;
static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u);
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w);

/**********************************************************************
    raise_g():

         -- calculates the contravariant form of a covariant tensor,
            using the inverse of the metric;
******************************************************************/
static void raise_g(CCTK_REAL vcov[NDIM], CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL vcon[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcon[i] = 0. ;
    for(j=0;j<NDIM;j++)
      vcon[i] += gcon[i][j]*vcov[j] ;
  }

  return ;
}

/**********************************************************************
     lower_g():

          -- calculates the ocvariant form of a contravariant tensor
             using the metric;
******************************************************************/
static void lower_g(CCTK_REAL vcon[NDIM], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL vcov[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcov[i] = 0. ;
    for(j=0;j<NDIM;j++)
      vcov[i] += gcov[i][j]*vcon[j] ;
  }

  return ;
}

/**********************************************************************
     ncov_calc():

         -- calculates the covariant form of the normal vector to our
            spacelike hypersurfaces ala the ADM formalism.

         -- requires the inverse metric;
******************************************************************/
static void ncov_calc(CCTK_REAL gcon[NDIM][NDIM],CCTK_REAL ncov[NDIM])
{
  CCTK_REAL lapse ;
  int i;

  lapse = sqrt(-1./gcon[0][0]) ;

  ncov[0] = -lapse ;
  for( i = 1; i < NDIM; i++) {
    ncov[i] = 0. ;
  }

  return ;
}

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/*
   pressure as a function of rho0 and u
   this is used by primtoU and Utoprim_?D
*/
static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u)
{
  DECLARE_CCTK_PARAMETERS;
  return((gamma_th /* <- Should be local polytropic Gamma factor */  - 1.)*u) ;
}



/*
   pressure as a function of rho0 and w = rho0 + u + p
   this is used by primtoU and Utoprim_1D
*/
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w)
{
  DECLARE_CCTK_PARAMETERS;
  return((gamma_th /* <- Should be local polytropic Gamma factor */ -1.)*(w - rho0)/gamma_th /* <- Should be local polytropic Gamma factor */ ) ;
}


