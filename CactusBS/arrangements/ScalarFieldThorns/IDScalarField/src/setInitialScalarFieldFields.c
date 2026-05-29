#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

/*
 * External functions
 */
void idscalarfield_setInitialScalarFields(CCTK_ARGUMENTS);


/*
 * Implementations
 */

void idscalarfield_setInitialScalarFields(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  /* get pointers to the scalar field variables (this avoids explicit
   * dependence on a scalar field thorn) */
  CCTK_REAL *scalarphi = CCTK_VarDataPtr(cctkGH, 0, scalarPhi);
  CCTK_REAL *scalarpi = CCTK_VarDataPtr(cctkGH, 0, scalarPi);
  if(NULL == scalarphi || NULL == scalarpi)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not find variables scalarPhi='%s' and scalarPi='%s'",
               scalarPhi, scalarPi);
    return; /* NOTREACHED*/
  }

  /* a scalar field is not conformally rescaled from the conformal initial data
   * space to the physical space */
#pragma omp parallel for
  for(int ind = 0 ; ind < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] ; ind++)
  {
    CCTK_REAL ini_phi, ini_grad_phi[3], ini_pi, ini_VV;

    calcInitialScalarFields(x[ind], y[ind], z[ind], &ini_phi, &ini_grad_phi[0], &ini_pi, &ini_VV);
    scalarphi[ind] = ini_phi;
    scalarpi[ind] = ini_pi;
  }
}
