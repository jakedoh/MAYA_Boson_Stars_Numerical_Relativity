#include "cctk.h"
#include "cctk_Arguments.h"

#ifdef __cplusplus
extern "C"
{
#endif
	void bondiid_getSourceInfo(double xx, double yy, double zz, double *rhohat, double jhat[3], double *psi, double *psiguess, double Aij[3][3]);
	double bondiid_getSourceMass(int index);
	double bondiid_getRestMass(int index);
	void bondiid_constructSources(CCTK_ARGUMENTS);
	CCTK_REAL bondiid_getExternalPsiGuess(CCTK_INT idx, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bondiid_getLorentzFactor(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bondiid_getConfDensity(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bondiid_getConfPress(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bondiid_getConfPolyConstant(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
#ifdef __cplusplus
}
#endif
