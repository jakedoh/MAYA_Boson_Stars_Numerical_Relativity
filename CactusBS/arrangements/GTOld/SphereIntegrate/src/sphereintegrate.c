#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

#include "sphereintegrate.h"

CCTK_INT SphereIntegrate_dtheta_dphi_setup(CCTK_INT ntheta, CCTK_INT nphi,
					   CCTK_REAL *dth, CCTK_REAL *dph,
					   CCTK_REAL *theta0, CCTK_REAL *phi0)
{
  DECLARE_CCTK_PARAMETERS;

  if (dth==NULL) {
    CCTK_WARN(1,"NULL pointer dth");
    return -1;
  }

  if (dph==NULL) {
    CCTK_WARN(1,"NULL pointer dph");
    return -1;
  }

  if (theta0==NULL) {
    CCTK_WARN(1,"NULL pointer theta0");
    return -1;
  }

  if (phi0==NULL) {
    CCTK_WARN(1,"NULL pointer phi0");
    return -1;
  }

  if (ntheta<2) {
    CCTK_WARN(1,"ntheta<2 - this won't work at all");
    return -1;
  }

  if (nphi<2) {
    CCTK_WARN(1,"nphi<2 - this won't work at all");
    return -1;
  }

  const CCTK_REAL pi=acos(-1.0);
  *dth=pi/(ntheta-1.);
  *dph=2.*pi/(nphi-1.);
  *theta0=0;
  *phi0=0;

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,
	       "coord setup: dth=%.19g dph=%.19g th0=%.19g ph0=%.19g",
	       *dth,*dph,*theta0,*phi0);
  }

  return 1;
}


CCTK_INT SphereIntegrate_setup_coords_open_theta(CCTK_INT ntheta, CCTK_INT nphi, CCTK_REAL rad,
						 CCTK_REAL origin_x, CCTK_REAL origin_y, CCTK_REAL origin_z,
						 CCTK_REAL *xx, CCTK_REAL *yy, CCTK_REAL *zz)
{
  DECLARE_CCTK_PARAMETERS;

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,
	       "SphereIntegrate_setup_coords_open_theta: r=%.19g origin=(%.19g,%.19g,%.19g)",
	       rad,origin_x,origin_y,origin_z);
  }

  if (xx==NULL) {
    CCTK_WARN(1,"NULL pointer in xx");
    return -1;
  }

  if (yy==NULL) {
    CCTK_WARN(1,"NULL pointer in yy");
    return -1;
  }

  if (zz==NULL) {
    CCTK_WARN(1,"NULL pointer in zz");
    return -1;
  }

  CCTK_REAL dth, dph, theta0, phi0;
  CCTK_REAL ierr=0;
  ierr=SphereIntegrate_dtheta_dphi_setup(ntheta, nphi, &dth, &dph, &theta0, &phi0);
  if (ierr<0) {
    CCTK_WARN(1,"failed to setup dtheta, dphi");
    return -1;
  }
  
  const CCTK_REAL dtp=dth*dph;

  const CCTK_INT npoints=ntheta*nphi;

  CCTK_REAL th,ct,st, ph,cp,sp;
  CCTK_INT ind2d;

  for (int i=0;i<npoints;i++) {
    xx[i]=0;yy[i]=0;zz[i]=0;
  }

  for (int i=0;i<ntheta;i++) {
    th=i*dth;
    ct=cos(th);
    st=sin(th);
    for (int j=0;j<nphi;j++) {
      ph=j*dph;
      cp=cos(ph);
      sp=sin(ph);
      ind2d=SPHEREINT_IND2D(i,j,nphi);
      xx[ind2d]=origin_x+ rad*cp*st;
      yy[ind2d]=origin_y+ rad*sp*st;
      zz[ind2d]=origin_z+ rad*   ct;
    }
  }

  return 1;
}


CCTK_INT SphereIntegrate_integrate_dOmega_open_theta(CCTK_REAL *integrand, 
						     CCTK_INT ntheta, CCTK_INT nphi,
						     CCTK_REAL *integral)
{
  DECLARE_CCTK_PARAMETERS;

  if (verbose>1) CCTK_INFO("in integrate_dOmega_open_theta");

  if (integrand==NULL) {
    CCTK_WARN(1,"NULL pointer in integrand for integrate_dOmega_open_theta!");
    return -1;
  }

  if (integral==NULL) {
    CCTK_WARN(1,"NULL pointer in integral for integrate_dOmega_open_theta!");
    return -1;
  }

  CCTK_REAL dth, dph, theta0, phi0;
  CCTK_REAL ierr=0;
  ierr=SphereIntegrate_dtheta_dphi_setup(ntheta, nphi, &dth, &dph, &theta0, &phi0);
  if (ierr<0) {
    CCTK_WARN(1,"failed to setup dtheta, dphi");
    return -1;
  }

  const CCTK_REAL dtp=dth*dph;

  CCTK_REAL th,sint,ph;
  CCTK_REAL iwtheta,iwphi,intweight,dOmega;
  CCTK_INT ind2d;

  CCTK_REAL intval,sum=0;
  CCTK_INT nan_count=0;

  for (int i=0;i<ntheta;i++) {
    th=i*dth;
    sint=sin(th);
    if      (i==0 || i==ntheta-1) iwtheta=0;
    else if (i==1 || i==ntheta-2) iwtheta=55.0/24.0;
    else if (i==2 || i==ntheta-3) iwtheta=-1.0/6.0;
    else if (i==3 || i==ntheta-4) iwtheta=11.0/8.0;
    else iwtheta=1;
    for (int j=0;j<nphi;j++) {
      ph=j*dph;
      if      (j==0 || j==nphi-1) iwphi=3.0/8.0;
      else if (j==1 || j==nphi-2) iwphi=7.0/6.0;
      else if (j==2 || j==nphi-3) iwphi=23.0/24.0;
      else iwphi=1;
      intweight=iwphi*iwtheta;

      dOmega=sint*dtp*intweight;

      ind2d=SPHEREINT_IND2D(i,j,nphi);

      intval=integrand[ind2d];
      if (i==0||i==ntheta-1) intval=0; // clean axis

      if (verbose>8) {
	CCTK_VInfo(CCTK_THORNSTRING,"i=%d j=%d ind2d=%d int=%.19g",
		   i,j,ind2d,intval);
      }
      
      // wipe NaN
      if (!isnan(intval)) {
	sum += intval*dOmega;
      }
      else {
	nan_count++;
	if (verbose>5) {
	  CCTK_VInfo(CCTK_THORNSTRING,"NaN at i=%d j=%d",i,j);
	}
      }
    }
  }
  *integral=sum;

  if (verbose>1) CCTK_VInfo(CCTK_THORNSTRING,"  integral=%.19g",sum);

  if (nan_count>0 && verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"found %d NaNs in integral",nan_count);
  }

  return 1;
}

