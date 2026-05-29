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
  CCTK_REAL *scalarphia = CCTK_VarDataPtr(cctkGH, 0, scalarPhia);
  CCTK_REAL *scalarpia = CCTK_VarDataPtr(cctkGH, 0, scalarPia);
  CCTK_REAL *scalarphib = CCTK_VarDataPtr(cctkGH, 0, scalarPhib);
  CCTK_REAL *scalarpib = CCTK_VarDataPtr(cctkGH, 0, scalarPib);

  if(NULL == scalarphia || NULL == scalarpia || NULL == scalarphib || NULL == scalarpib)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not find variables scalarPhia='%s', scalarPia='%s', scalarPhib='%s', and scalarPib='%s'",
               scalarPhia, scalarPia,scalarPhib, scalarPib);
    return; /* NOTREACHED*/
  }

  /* a scalar field is not conformally rescaled from the conformal initial data
   * space to the physical space */
#pragma omp parallel for
  for(int ind = 0 ; ind < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] ; ind++)
  { 
    CCTK_REAL ini_phia, ini_grad_phia[3], ini_pia, ini_phib, ini_grad_phib[3], ini_pib, ini_gvphiaa, ini_gvphibb, ini_gvphiab, ini_VV;
    calcInitialScalarFields(x[ind], y[ind], z[ind], &ini_phia, &ini_grad_phia[0], &ini_pia, &ini_phib, &ini_grad_phib[0], &ini_pib, &ini_VV);
    scalarphia[ind] = ini_phia;
    scalarpia[ind] = ini_pia;
    scalarphib[ind] = ini_phib;
    scalarpib[ind] = ini_pib;
  }
}
