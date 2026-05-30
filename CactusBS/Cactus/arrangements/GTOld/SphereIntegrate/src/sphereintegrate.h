#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 


/* access to the points on the sphere is by convention */
#define SPHEREINT_IND2D(i, j, nphi)    ( (j) + (i)*(nphi) )


CCTK_INT SphereIntegrate_dtheta_dphi_setup(CCTK_INT ntheta, CCTK_INT nphi,
					   CCTK_REAL *dth, CCTK_REAL *dph,
					   CCTK_REAL *theta0, CCTK_REAL *phi0);

CCTK_INT SphereIntegrate_setup_coords_open_theta(CCTK_INT ntheta, CCTK_INT nphi, CCTK_REAL rad,
						 CCTK_REAL origin_x, CCTK_REAL origin_y, CCTK_REAL origin_z,
						 CCTK_REAL *xx, CCTK_REAL *yy, CCTK_REAL *zz);

CCTK_INT SphereIntegrate_integrate_dOmega_open_theta(CCTK_REAL *integrand, 
						     CCTK_INT ntheta, CCTK_INT nphi,
						     CCTK_REAL *integral);

