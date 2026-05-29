#include <stdio.h>
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

void fillSymmetric4by4Matrix(double mat[4][4])
{
	mat[1][0] = mat[0][1];
	mat[2][0] = mat[0][2];
	mat[3][0] = mat[0][3];
	mat[2][1] = mat[1][2];
	mat[3][1] = mat[1][3];
	mat[3][2] = mat[2][3];
}

double calcJi(int idx, double normal[4], double energy_momentum[4][4], double inverse_metric[4][4])
{
	double ji = 0.;
	for(int aa = 0; aa < 4; aa++)
	{
		for(int kk = 1; kk < 4; kk++)
		{
		ji += normal[aa]*energy_momentum[kk][aa]*inverse_metric[kk][idx];
		}
	}
	return ji;
}

void MatterAnalysis_CalculateQuantities(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_LOOP3_ALL( calc_matteranalysis_quantities, cctkGH, i, j, k)
  {
        CCTK_INT ind=CCTK_GFINDEX3D(cctkGH,i,j,k);
	double energy_momentum[4][4];
	energy_momentum[0][0] = eTtt[ind];
	energy_momentum[0][1] = eTtx[ind];
	energy_momentum[0][2] = eTty[ind];
	energy_momentum[0][3] = eTtz[ind];
	energy_momentum[1][1] = eTxx[ind];
	energy_momentum[1][2] = eTxy[ind];
	energy_momentum[1][3] = eTxz[ind];
	energy_momentum[2][2] = eTyy[ind];
	energy_momentum[2][3] = eTyz[ind];
	energy_momentum[3][3] = eTzz[ind];

	fillSymmetric4by4Matrix(energy_momentum);

	double inverse_metric[4][4];

	double lapse, shiftx, shifty, shiftz;
	lapse = alp[ind];
	shiftx = betax[ind];
	shifty = betay[ind];
	shiftz = betaz[ind];

	inverse_metric[0][0] = -1./(lapse*lapse);
	inverse_metric[0][1] = shiftx/(lapse*lapse);
	inverse_metric[0][2] = shifty/(lapse*lapse);
	inverse_metric[0][3] = shiftz/(lapse*lapse);

	double metric[4][4];
	metric[0][0] = -lapse*lapse + shiftx*shiftx + shifty*shifty + shiftz*shiftz;
	metric[0][1] = shiftx;
	metric[0][2] = shifty;
	metric[0][3] = shiftz;
	metric[1][1] = gxx[ind];
	metric[1][2] = gxy[ind];
	metric[1][3] = gxz[ind];
	metric[2][2] = gyy[ind];
	metric[2][3] = gyz[ind];
	metric[3][3] = gzz[ind];

	fillSymmetric4by4Matrix(metric);

	double spatial_determinant = 2*metric[1][2]*metric[1][3]*metric[2][3] + metric[3][3]*(metric[1][1]*metric[2][2] - metric[2][1]*metric[2][1]) - metric[2][2]*(metric[3][1]*metric[3][1]) - metric[1][1]*metric[3][2]*metric[3][2];

	inverse_metric[1][1] = (metric[2][2]*metric[3][3] - metric[2][3]*metric[3][2])/spatial_determinant;
	inverse_metric[1][2] = (metric[1][3]*metric[3][2] - metric[1][2]*metric[3][3])/spatial_determinant;
	inverse_metric[1][3] = (metric[1][2]*metric[2][3] - metric[1][3]*metric[2][2])/spatial_determinant;
	inverse_metric[2][2] = (metric[1][1]*metric[3][3] - metric[1][3]*metric[3][1])/spatial_determinant;
	inverse_metric[2][3] = (metric[1][3]*metric[2][1] - metric[1][1]*metric[2][3])/spatial_determinant;
	inverse_metric[3][3] = (metric[1][1]*metric[2][2] - metric[1][2]*metric[2][1])/spatial_determinant;

	fillSymmetric4by4Matrix(inverse_metric);

	double normal[4];
	normal[0] = 1./lapse;
	normal[1] = shiftx/lapse;
	normal[2] = shifty/lapse;
	normal[3] = shiftz/lapse;

	//density_weighted_position stuff
	rho_times_x[ind] = x[ind]*rho[ind];
	rho_times_y[ind] = y[ind]*rho[ind];
	rho_times_z[ind] = z[ind]*rho[ind];

	double conformal_factor = pow(spatial_determinant, 1./12.);

	eJx[ind] = calcJi(1, normal, energy_momentum, inverse_metric)*pow(conformal_factor, 4.0);
	eJy[ind] = calcJi(2, normal, energy_momentum, inverse_metric)*pow(conformal_factor, 4.0);
	eJz[ind] = calcJi(3, normal, energy_momentum, inverse_metric)*pow(conformal_factor, 4.0);
  } CCTK_ENDLOOP3_ALL( calc_matteranalysis_quantities );
}

