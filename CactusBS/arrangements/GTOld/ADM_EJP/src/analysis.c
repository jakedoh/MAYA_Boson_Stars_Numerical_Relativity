#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

#include "sphereintegrate.h"

static int ADM_EJP_write_data(CCTK_ARGUMENTS,
                              CCTK_INT det,
                              CCTK_REAL rad);

static CCTK_INT have_integrand_memory=0;
static CCTK_REAL *Eintegrand;
static CCTK_REAL *Pxintegrand;  
static CCTK_REAL *Pyintegrand;  
static CCTK_REAL *Pzintegrand;  
static CCTK_REAL *Jxintegrand;  
static CCTK_REAL *Jyintegrand;  
static CCTK_REAL *Jzintegrand;  
static CCTK_REAL *dEdtintegrand;
static CCTK_REAL *dPxdtintegrand;
static CCTK_REAL *dPydtintegrand;
static CCTK_REAL *dPzdtintegrand;
static CCTK_REAL *dSxdtintegrand;
static CCTK_REAL *dSydtintegrand;
static CCTK_REAL *dSzdtintegrand;

static int allocate_memory(CCTK_INT npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>1) CCTK_INFO("in allocate_memory");
    

  if (!have_integrand_memory) {
    if (verbose>1) CCTK_INFO("allocating new memory");
    Eintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Eintegrand==NULL) {
      CCTK_WARN(1,"Eintegrand allocation failed");
      return -1;
    }
    Pxintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Pxintegrand==NULL) {
      CCTK_WARN(1,"Pxintegrand allocation failed");
      return -1;
    }
    Pyintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Pyintegrand==NULL) {
      CCTK_WARN(1,"Pyintegrand allocation failed");
      return -1;
    }
    Pzintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Pzintegrand==NULL) {
      CCTK_WARN(1,"Pzintegrand allocation failed");
      return -1;
    }
    Jxintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Jxintegrand==NULL) {
      CCTK_WARN(1,"Jxintegrand allocation failed");
      return -1;
    }
    Jyintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Jyintegrand==NULL) {
      CCTK_WARN(1,"Jyintegrand allocation failed");
      return -1;
    }
    Jzintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (Jzintegrand==NULL) {
      CCTK_WARN(1,"Jzintegrand allocation failed");
      return -1;
    }

    // make sure nothing touches the memory later if not requested!
    if (compute_gw_rad) {

      dEdtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dEdtintegrand==NULL) {
	CCTK_WARN(1,"dEdtintegrand allocation failed");
	return -1;
      }
      dPxdtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dPxdtintegrand==NULL) {
	CCTK_WARN(1,"dPxdtintegrand allocation failed");
	return -1;
      }
      dPydtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dPydtintegrand==NULL) {
	CCTK_WARN(1,"dPydtintegrand allocation failed");
	return -1;
      }
      dPzdtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dPzdtintegrand==NULL) {
	CCTK_WARN(1,"dPzdtintegrand allocation failed");
	return -1;
      }
      dSxdtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dSxdtintegrand==NULL) {
	CCTK_WARN(1,"dSxdtintegrand allocation failed");
	return -1;
      }
      dSydtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dSydtintegrand==NULL) {
	CCTK_WARN(1,"dSydtintegrand allocation failed");
	return -1;
      }
      dSzdtintegrand = (CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
      if (dSzdtintegrand==NULL) {
	CCTK_WARN(1,"dSzdtintegrand allocation failed");
	return -1;
      }
    }

    have_integrand_memory=1;
    return 1;
  }
  else {
    if (verbose>1) CCTK_INFO("already allocated memory");
    return 2;
  }
  return 1;
}

void ADM_EJP_Compute_EJP (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ind=0;
  CCTK_REAL ierr=0;


  if (verbose>1) {
    CCTK_INFO("in compute_ejp");
  }

  if (*do_nothing==1) {
    return;
  }

  if (*rdet<=0) {
    // update detector
    *current_detector=*current_detector-1;
    CCTK_VInfo(CCTK_THORNSTRING,
               "det rad negative in analysis detector: %d number_of_dets: %d rdet=%g",
               *current_detector,number_of_detectors,*rdet);
    return;
  }

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"analysis.c: detector no. %d",*current_detector);
  }

  const CCTK_INT det=*current_detector;
  const CCTK_REAL rad=*rdet;

  // 2D surface integral
  CCTK_REAL th=0,ph=0,rP4=0,iP4=0,absP4sq=0;
  CCTK_REAL sint=0,cost=0,sinp=0,cosp=0;
  CCTK_REAL nx=0,ny=0,nz=0;
  CCTK_REAL dP4rdphi,dP4idphi;

  // surface integrals
  CCTK_REAL Eint=0;
  CCTK_REAL Pxint=0,Pyint=0,Pzint=0;
  CCTK_REAL Jxint=0,Jyint=0,Jzint=0;
  CCTK_INT ind2d,indp,indm;
  CCTK_REAL ux,uy,uz;
  CCTK_REAL term1,term2,term3;

  // adm mass quantities
  CCTK_REAL dxdxx,dxdxy,dxdxz,dxdyy,dxdyz,dxdzz, dxdyx,dxdzx,dxdzy;
  CCTK_REAL dydxx,dydxy,dydxz,dydyy,dydyz,dydzz, dydyx,dydzx,dydzy;
  CCTK_REAL dzdxx,dzdxy,dzdxz,dzdyy,dzdyz,dzdzz, dzdyx,dzdzx,dzdzy;

  // adm linear momentum quantities
  CCTK_REAL dxx,dxy,dxz,dyy,dyz,dzz, dyx,dzx,dzy;
  CCTK_REAL exx,exy,exz,eyy,eyz,ezz, eyx,ezx,ezy;
  CCTK_REAL uxx,uxy,uxz,uyy,uyz,uzz, uyx,uzx,uzy;
  CCTK_REAL detm,idetm;
  CCTK_REAL kdiagx,kdiagy,kdiagz,skdiag;
  CCTK_REAL nlx,nly,nlz;

  // adm angular momentum quantities
  CCTK_REAL kxa_na,kya_na,kza_na;

  // T00
  CCTK_REAL bx,by,bz, dxbx,dybx,dzbx, dxby,dyby,dzby, dxbz,dybz,dzbz;
  CCTK_REAL lbgxx,lbgxy,lbgxz,lbgyy,lbgyz,lbgzz;
  CCTK_REAL pf, dtgxx,dtgxy,dtgxz,dtgyx,dtgyy,dtgyz,dtgzx,dtgzy,dtgzz;
  CCTK_REAL dEdtdOmega=0, dEdt=0,dPxdt=0,dPydt=0,dPzdt=0,dSxdt=0,dSydt=0,dSzdt=0;

  CCTK_REAL iwtheta, iwphi, intweight;

  // setup integrands
  const CCTK_INT npoints=ntheta*nphi;

  ierr=allocate_memory(npoints);
  if (ierr<0) {
    CCTK_WARN(1,"memory allocation failed");
    return;
  }
  if(ierr==2 && verbose>1) {
    CCTK_INFO("memory already allocated");
  }


  CCTK_REAL dth, dph, theta0, phi0;
  ierr=SphereIntegrate_dtheta_dphi_setup(ntheta, nphi, &dth, &dph, &theta0, &phi0);
  if (ierr<0) {
    CCTK_WARN(1,"failed to setup dtheta, dphi");
    return;
  }

  CCTK_REAL rad2=rad*rad;

  for (int i=0;i<ntheta;i++) {
    th=theta0+i*dth;
    sint=sin(th);
    cost=cos(th);

    for (int j=0;j<nphi;j++) {
      ph=phi0+j*dph;
      sinp=sin(ph);
      cosp=cos(ph);

      ind2d=j+nphi*i;

      // normal vector
      ux=cosp*sint;
      uy=sinp*sint;
      uz=     cost;

      dxdxx=dx_gxx[ind2d];
      dxdxy=dx_gxy[ind2d];
      dxdxz=dx_gxz[ind2d];
      dxdyy=dx_gyy[ind2d];
      dxdyz=dx_gyz[ind2d];
      dxdzz=dx_gzz[ind2d];
      dxdyx=dxdxy;
      dxdzx=dxdxz;
      dxdzy=dxdyz;

      dydxx=dy_gxx[ind2d];
      dydxy=dy_gxy[ind2d];
      dydxz=dy_gxz[ind2d];
      dydyy=dy_gyy[ind2d];
      dydyz=dy_gyz[ind2d];
      dydzz=dy_gzz[ind2d];
      dydyx=dydxy;
      dydzx=dydxz;
      dydzy=dydyz;

      dzdxx=dz_gxx[ind2d];
      dzdxy=dz_gxy[ind2d];
      dzdxz=dz_gxz[ind2d];
      dzdyy=dz_gyy[ind2d];
      dzdyz=dz_gyz[ind2d];
      dzdzz=dz_gzz[ind2d];
      dzdyx=dzdxy;
      dzdzx=dzdxz;
      dzdzy=dzdyz;

      term1=(dxdxx+dydyx+dzdzx)-(dxdxx+dxdyy+dxdzz);
      term2=(dxdxy+dydyy+dzdzy)-(dydxx+dydyy+dydzz);
      term3=(dxdxz+dydyz+dzdzz)-(dzdxx+dzdyy+dzdzz);

      // energy
      // integral over dA=r^2*dOmega;
      Eintegrand[ind2d] = (ux*term1+uy*term2+uz*term3)*rad2;
      if (i==0||i==ntheta-1) Eintegrand[ind2d] = 0;

      // down-stairs metric
      dxx=gxx2d[ind2d];
      dxy=gxy2d[ind2d];
      dxz=gxz2d[ind2d];
      dyy=gyy2d[ind2d];
      dyz=gyz2d[ind2d];
      dzz=gzz2d[ind2d];
      dyx=dxy;
      dzx=dxz;
      dzy=dyz;


      // determinant of 3-metric
      detm=dxx*dyy*dzz + 2.*dxy*dxz*dyz
        -(dxx*dyz*dyz+dyy*dxz*dxz+dzz*dxy*dxy);
      idetm=1./detm;

      // up-stairs metric
      uxx=idetm*(dyy*dzz-dyz*dyz);
      uxy=idetm*(dxz*dyz-dzz*dxy);
      uxz=idetm*(dxy*dyz-dyy*dxz);
      uyy=idetm*(dxx*dzz-dxz*dxz);
      uyz=idetm*(dxy*dxz-dxx*dyz);
      uzz=idetm*(dxx*dyy-dxy*dxy);
      uyx=uxy;
      uzx=uxz;
      uzy=uyz;

      // extr. curv
      exx=kxx2d[ind2d];
      exy=kxy2d[ind2d];
      exz=kxz2d[ind2d];
      eyy=kyy2d[ind2d];
      eyz=kyz2d[ind2d];
      ezz=kzz2d[ind2d];
      eyx=exy;
      ezx=exz;
      ezy=eyz;

      // K^x_x, K^y_y, K^z_z
      kdiagx=uxx*exx+uxy*eyx+uxz*ezx;
      kdiagy=uyx*exy+uyy*eyy+uyz*ezy;
      kdiagz=uzx*exz+uzy*eyz+uzz*ezz;
      skdiag=kdiagx+kdiagy+kdiagz;

      // linear momentum
      // down-stairs unit normal
      nlx=dxx*ux+dxy*uy+dxz*uz;
      nly=dyx*ux+dyy*uy+dyz*uz;
      nlz=dzx*ux+dzy*uy+dzz*uz;
      Pxintegrand[ind2d] = (exx*ux+eyx*uy+ezx*uz +skdiag*nlx)*rad2;
      Pyintegrand[ind2d] = (exy*ux+eyy*uy+ezy*uz +skdiag*nly)*rad2;
      Pzintegrand[ind2d] = (exz*ux+eyz*uy+ezz*uz +skdiag*nlz)*rad2;
      if (i==0||i==ntheta-1) {
	Pxintegrand[ind2d] = 0;
	Pyintegrand[ind2d] = 0;
	Pzintegrand[ind2d] = 0;
      }

      // angular momentum
      kxa_na = ux*exx+uy*exy+uz*exz;
      kya_na = ux*eyx+uy*eyy+uz*eyz;
      kza_na = ux*ezx+uy*ezy+uz*ezz;
      term1 = uzx*kxa_na
             +uzy*kya_na
             +uzz*kza_na;
      term2 = uyx*kxa_na
             +uyy*kya_na
             +uyz*kza_na;
      Jxintegrand[ind2d] = (uy*term1-uz*term2)*rad*rad2;
      if (i==0||i==ntheta-1) Jxintegrand[ind2d] = 0;


      term1 = uxx*kxa_na
             +uxy*kya_na
             +uxz*kza_na;
      term2 = uzx*kxa_na
             +uzy*kya_na
             +uzz*kza_na;
      Jyintegrand[ind2d] = (uz*term1-ux*term2)*rad*rad2;
      if (i==0||i==ntheta-1) Jyintegrand[ind2d] = 0;

      term1 = uyx*kxa_na
             +uyy*kya_na
             +uyz*kza_na;
      term2 = uxx*kxa_na
             +uxy*kya_na
             +uxz*kza_na;
      Jzintegrand[ind2d] = (ux*term1-uy*term2)*rad*rad2;
      if (i==0||i==ntheta-1) Jzintegrand[ind2d] = 0;

      // GW T_00
      if (compute_gw_rad)
      {
	bx  =   betax2d[ind2d];
	dxbx=dx_betax2d[ind2d];
	dybx=dy_betax2d[ind2d];
	dzbx=dz_betax2d[ind2d];
	by  =   betay2d[ind2d];
	dxby=dx_betay2d[ind2d];
	dyby=dy_betay2d[ind2d];
	dzby=dz_betay2d[ind2d];
	bz  =   betaz2d[ind2d];
	dxbz=dx_betaz2d[ind2d];
	dybz=dy_betaz2d[ind2d];
	dzbz=dz_betaz2d[ind2d];
	
	lbgxx = bx*dzdxx+bx*dydxx+2*dxbz*dxz+2*dxby*dxy+2*dxbx*dxx+bx*dxdxx;
	lbgxy = bx*dzdxy+dxbz*dyz+dxby*dyy+bx*dydxy+dxz*dybz+dxy*dyby+dxx 
	  *dybx+dxbx*dxy+bx*dxdxy;
	lbgxz = dxbz*dzz+bx*dzdxz+dxz*dzbz+dxy*dzby+dxx*dzbx+dxby*dyz+bx*dydxz 
	  +dxbx*dxz+bx*dxdxz;
	lbgyy = by*dzdyy+2*dybz*dyz+2*dyby*dyy+by*dydyy+2*dxy*dybx+by*dxdyy;
	lbgyz = dybz*dzz+by*dzdyz+dyz*dzbz+dyy*dzby+dxy*dzbx+dyby*dyz+by*dydyz
	  +dxz*dybx+by*dxdyz;
	lbgzz = 2*dzbz*dzz+bz*dzdzz+2*dyz*dzby+2*dxz*dzbx+bz*dydzz+bz*dxdzz;
	
	pf = -2.*alpha2d[ind2d];
	dtgxx = pf*exx+lbgxx;
	dtgxy = pf*exy+lbgxy;
	dtgxz = pf*exz+lbgxz;
	dtgyy = pf*eyy+lbgyy;
	dtgyz = pf*eyz+lbgyz;
	dtgzz = pf*ezz+lbgzz;
	dtgyx = dtgxy;
	dtgzx = dtgxz;
	dtgzy = dtgyz;
	

	term1=dtgxx+2.*dtgxy+2.*dtgxz+dtgyy+2.*dtgyz+dtgzz;
	dEdtdOmega = term1*term1;
	
	dEdtintegrand[ind2d] = dEdtdOmega*rad2;
	if (i==0||i==ntheta-1) dEdtintegrand[ind2d] = 0;

	dPxdtintegrand[ind2d] = nlx*dEdtdOmega*rad2;
	dPydtintegrand[ind2d] = nly*dEdtdOmega*rad2;
	dPzdtintegrand[ind2d] = nlz*dEdtdOmega*rad2;
	if (i==0||i==ntheta-1) {
	  dPxdtintegrand[ind2d] = 0;
	  dPydtintegrand[ind2d] = 0;
	  dPzdtintegrand[ind2d] = 0;
	}

	term1 = 0.5*((dtgzz*dydzz+2*dtgyz*dydyz+dtgyy*dydyy+2*dtgxz*dydxz
		      +2*dtgxy*dydxy+dtgxx*dydxx)*nlz+(-dtgzz*dzdzz-2*dtgyz*dzdyz-dtgyy
						       *dzdyy-2*dtgxz*dzdxz-2*dtgxy*dzdxy-dtgxx*dzdxx)*nly-2*dtgyz*dzz
		     +(2*dtgzz-2*dtgyy)*dyz+2*dtgyz*dyy-2*dtgxy*dxz+2*dtgxz*dxy);

	dSxdtintegrand[ind2d] = term1*rad2;
	if (i==0||i==ntheta-1) dSxdtintegrand[ind2d] = 0;
	
	term1 = -0.5*((dtgzz*dxdzz+2*dtgyz*dxdyz+dtgyy*dxdyy+2*dtgxz*dxdxz
		       +2*dtgxy*dxdxy+dtgxx*dxdxx)*nlz+(-dtgzz*dzdzz-2*dtgyz*dzdyz-dtgyy
							*dzdyy-2*dtgxz*dzdxz-2*dtgxy*dzdxy-dtgxx*dzdxx)*nlx-2*dtgxz*dzz
		      -2*dtgxy*dyz+(2*dtgzz-2*dtgxx)*dxz+2*dtgyz*dxy+2*dtgxz*dxx);
	
	dSydtintegrand[ind2d] = term1*rad2;
	if (i==0||i==ntheta-1) dSydtintegrand[ind2d] = 0;
	
	term1 = 0.5*((dtgzz*dxdzz+2*dtgyz*dxdyz+dtgyy*dxdyy+2*dtgxz*dxdxz
		      +2*dtgxy*dxdxy+dtgxx*dxdxx)*nly+(-dtgzz*dydzz-2*dtgyz*dydyz-dtgyy
						       *dydyy-2*dtgxz*dydxz-2*dtgxy*dydxy-dtgxx*dydxx)*nlx-2*dtgxz*dyz
		     -2*dtgxy*dyy+2*dtgyz*dxz+(2*dtgyy-2*dtgxx)*dxy+2*dtgxy*dxx);

	dSzdtintegrand[ind2d] = term1*rad2;
	if (i==0||i==ntheta-1) dSzdtintegrand[ind2d] = 0;
      }
    }
  }

  // do integrals over dOmega
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Eintegrand, ntheta, nphi, &Eint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do E integral");
  }
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Pxintegrand, ntheta, nphi, &Pxint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do Px integral");
  }
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Pyintegrand, ntheta, nphi, &Pyint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do Py integral");
  }
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Pzintegrand, ntheta, nphi, &Pzint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do Pz integral");
  }
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Jxintegrand, ntheta, nphi, &Jxint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do Jx integral");
  }
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Jyintegrand, ntheta, nphi, &Jyint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do Jy integral");
  }
  ierr=SphereIntegrate_integrate_dOmega_open_theta(Jzintegrand, ntheta, nphi, &Jzint);
  if (ierr<0) {
    CCTK_WARN(1,"failed to do Jz integral");
  }

  // fix normalization
  const CCTK_REAL pi=acos(-1.0);
  const CCTK_REAL prefac=1./pi;

  // E(t)
  energy[det]=prefac/16. * Eint;

  // P^i(t)
  Px[det]=prefac/8. * Pxint;
  Py[det]=prefac/8. * Pyint;
  Pz[det]=prefac/8. * Pzint;

  // J_i(t)
  Jx[det]=prefac/8. * Jxint;
  Jy[det]=prefac/8. * Jyint;
  Jz[det]=prefac/8. * Jzint;

  if (compute_gw_rad) {
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dEdtintegrand, ntheta, nphi, &dEdt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dEdt integral");
    }
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dPxdtintegrand, ntheta, nphi, &dPxdt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dPxdt integral");
    }
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dPydtintegrand, ntheta, nphi, &dPydt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dPydt integral");
    }
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dPzdtintegrand, ntheta, nphi, &dPzdt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dPzdt integral");
    }
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dSxdtintegrand, ntheta, nphi, &dSxdt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dSxdt integral");
    }
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dSydtintegrand, ntheta, nphi, &dSydt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dSydt integral");
    }
    ierr=SphereIntegrate_integrate_dOmega_open_theta(dSzdtintegrand, ntheta, nphi, &dSzdt);
    if (ierr<0) {
      CCTK_WARN(1,"failed to do dSzdt integral");
    }
  
    // E_rad
    const CCTK_REAL deltat=cctk_time- lastt[*current_detector];
    Erad[det]+=prefac/32.*deltat*dEdt;

    // P_rad
    Pxrad[det]+=prefac/32.*deltat*dPxdt;
    Pyrad[det]+=prefac/32.*deltat*dPydt;
    Pzrad[det]+=prefac/32.*deltat*dPzdt;
    
    // S_rad
    Sxrad[det]+=prefac/16.*deltat*dSxdt;
    Syrad[det]+=prefac/16.*deltat*dSydt;
    Szrad[det]+=prefac/16.*deltat*dSzdt;
  }
  else {
    Erad[det]=0;
    Pxrad[det]=0;
    Pyrad[det]=0;
    Pzrad[det]=0;
    Sxrad[det]=0;
    Syrad[det]=0;
    Szrad[det]=0;
  }

  // IO
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"r=%g E=%g P^i=(%g,%g,%g) J=(%g,%g,%g)",rad,
               energy[det], Px[det],Py[det],Pz[det], Jx[det],Jy[det],Jz[det]);
    if(compute_gw_rad) {
      CCTK_VInfo(CCTK_THORNSTRING,"  E_rad=%g P_rad^i=(%g,%g,%g) J_rad=(%g,%g,%g)",
		 Erad[det], Pxrad[det],Pyrad[det],Pzrad[det], Sxrad[det],Syrad[det],Szrad[det]);
    }
  }

  // write data only on CPU 0
  const CCTK_INT myproc=CCTK_MyProc(cctkGH);
  if (myproc==0) {
    CCTK_INT ierr=ADM_EJP_write_data(CCTK_PASS_CTOC,
                                     det,
                                     rad);
    if (ierr<0) {
      CCTK_WARN(1,"data writing failed");
    }
  }
}



static int ADM_EJP_write_data(CCTK_ARGUMENTS,
                              CCTK_INT det,
                              CCTK_REAL rad)                              
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char fmode[1024];
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  // file mode: append if already written
  if (*adm_ejp_files_exist>0) {
    sprintf(fmode,"a");
  }
  else {
    sprintf(fmode,"w");
  }

  // filename
  sprintf(varname, "adm_ejp_det_%d",det);

  filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                              strlen (file_extension) +2);
  sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "adm_ejp_write: Could not open scalar output file '%s'",
                filename);
    return -1;
  }

  // write header on startup
  if (*adm_ejp_files_exist<=0) {
    fprintf(file,"# ADM EJP Quantities\n");
    fprintf(file,"# det no.=%d rdet=%g\n",det,rad);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:t 2:E 3:Px 4:Py 5:Pz 6:Jx 7:Jy 8:Jz\n");
  }

  // write data
  sprintf (format_str_real,
           "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
           out_format,out_format,out_format,out_format,out_format,out_format,
           out_format,out_format);
  fprintf(file, format_str_real,
	  cctk_time, energy[det],
          Px[det], Py[det], Pz[det],
          Jx[det], Jy[det], Jz[det]);

  fclose(file); // end of ADM files
  free(filename); // free filename, so that the pointer can be used below

  // write an extra file for TT radiation quantities
  if(compute_gw_rad) {
    // filename - allocate again, because length is different
    sprintf(varname, "tt_rad_det_%d",det);

    filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
				strlen (file_extension) +2);
    sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);
    
    // open file
    file = fopen (filename, fmode);
    if (!file) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "adm_ejp_write: Could not open scalar output file '%s'",
		  filename);
      return -1;
    }
    
    // write header on startup
    if (*adm_ejp_files_exist<=0) {
      fprintf(file,"# TT radiation Quantities\n");
      fprintf(file,"# det no.=%d rdet=%g\n",det,rad);
      fprintf(file,"# gnuplot column index:\n");
      fprintf(file,"# 1:t 2:Erad 3:Pxrad 4:Pyrad 5:Pzrad 6:Jxrad 7:Jyrad 8:Jzrad\n");
    }
    
    // write data
    sprintf (format_str_real,
	     "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
	     out_format,out_format,out_format,out_format,out_format,out_format,
	     out_format,out_format);
    fprintf(file, format_str_real, 
	    cctk_time, Erad[det],
	    Pxrad[det],Pyrad[det],Pzrad[det],
	    Sxrad[det],Syrad[det],Szrad[det]);
    fclose(file);
    free(filename); // free filename again
  }


  if (det==number_of_detectors-1) { // we have written all detectors
    *adm_ejp_files_exist=1;
  }
  return 1;
}
