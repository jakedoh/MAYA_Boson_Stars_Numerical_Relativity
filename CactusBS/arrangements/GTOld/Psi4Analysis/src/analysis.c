#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 


static int psi4analysis_write_data(CCTK_ARGUMENTS,
                                   CCTK_INT det,
                                   CCTK_REAL rad,
                                   CCTK_REAL dEdt,
                                   CCTK_REAL dPxdt, CCTK_REAL dPydt, CCTK_REAL dPzdt,
				   CCTK_REAL dJxdt, CCTK_REAL dJydt, CCTK_REAL dJzdt);


void Psi4Analysis_Compute_EPRad (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ind=0;

  if (verbose>1) {
    CCTK_INFO("in compute_eprad");
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

  // delta t
  const CCTK_REAL deltat=cctk_time - lastt[*current_detector];

  const CCTK_REAL rad=*rdet;
  
  // 2D surface integral
  CCTK_REAL dOmega=0;
  const CCTK_REAL pi=acos(-1.0);
  const CCTK_REAL dth=pi/(*ntheta-1);
  const CCTK_REAL dph=2.*pi/(*nphi-1);
  const CCTK_REAL idph=1./dph;
  const CCTK_REAL idth=1./dth;
  const CCTK_REAL i2dph=idph/2.;
  const CCTK_REAL i2dth=idth/2.;
  const CCTK_REAL dtp=dth*dph;
  
  CCTK_REAL th=0,ph=0,rP4=0,iP4=0,absP4sq=0;
  CCTK_REAL sint=0,cost=0,sinp=0,cosp=0,cott=0,csct=0;
  CCTK_REAL nx=0,ny=0,nz=0;
  CCTK_REAL xv=0,yv=0,zv=0;
  CCTK_REAL dP4rdphi,dP4idphi;
  CCTK_REAL iwtheta=0,iwphi=0,intweight=0;

  CCTK_REAL rdphIP4dt, idphIP4dt;
  CCTK_REAL rdthIP4dt, idthIP4dt;

  // surface integrals
  CCTK_REAL dEdt=0;
  CCTK_REAL dPxdt=0,dPydt=0,dPzdt=0;
  CCTK_REAL dPrdp,dPidp,dPrdt,dPidt,I2Pr,I2Pi;
  CCTK_REAL dJxdt=0,dJydt=0,dJzdt=0;
  CCTK_REAL f1=0,f2=0;
  CCTK_INT ind2d,ind3d,indp,indm;
  for (int i=0;i<*ntheta;i++) {
    th=i*dth;
    sint=sin(th);
    cost=cos(th);
    if (i>0) {
      cott=1/tan(th);
      csct=1/sint;
    }
    else {
      cott=0;
      csct=0;
    }

    if      (i==0 || i==*ntheta-1) iwtheta=0;
    else if (i==1 || i==*ntheta-2) iwtheta=55.0/24.0;
    else if (i==2 || i==*ntheta-3) iwtheta=-1.0/6.0;
    else if (i==3 || i==*ntheta-4) iwtheta=11.0/8.0;
    else iwtheta=1;
    for (int j=0;j<*nphi;j++) {
      ph=j*dph;
      sinp=sin(ph);
      cosp=cos(ph);
      if      (j==0 || j==*nphi-1) iwphi=3.0/8.0;
      else if (j==1 || j==*nphi-2) iwphi=7.0/6.0;
      else if (j==2 || j==*nphi-3) iwphi=23.0/24.0;
      else iwphi=1;
      intweight=iwphi*iwtheta;

      dOmega=sint*dtp*intweight;

      ind3d=i+ maxntheta*j+ maxntheta* maxnphi*(*current_detector);

      rP4=int_psi4r_dt[ind3d];
      iP4=int_psi4i_dt[ind3d];

      absP4sq=rP4*rP4+iP4*iP4;

      // radiated quantities
      // energy
      dEdt+=absP4sq*dOmega;

      nx=cosp*sint;
      ny=sinp*sint;
      nz=     cost;

      // linear momentum
      dPxdt+=nx*absP4sq*dOmega;
      dPydt+=ny*absP4sq*dOmega;
      dPzdt+=nz*absP4sq*dOmega;

      // angular momentum
      // d_ph I(P4)dt
      if (j>0 && j<*nphi-1) {
	indp=i+ maxntheta*(j+1)+ maxntheta* maxnphi*(*current_detector);
	indm=i+ maxntheta*(j-1)+ maxntheta* maxnphi*(*current_detector);
	rdphIP4dt=(int_psi4r_dt[indp]-int_psi4r_dt[indm])*i2dph;
	idphIP4dt=(int_psi4i_dt[indp]-int_psi4i_dt[indm])*i2dph;
      }
      else if (j==0) {
	indp=i+ maxntheta*(1)+ maxntheta* maxnphi*(*current_detector); 
	indm=i+ maxntheta*(0)+ maxntheta* maxnphi*(*current_detector);
	rdphIP4dt=(int_psi4r_dt[indp]-int_psi4r_dt[indm])*idph;
	idphIP4dt=(int_psi4i_dt[indp]-int_psi4i_dt[indm])*idph;
      }
      else if (j==*nphi-1) {
	indp=i+ maxntheta*(*nphi-1)+ maxntheta* maxnphi*(*current_detector); 
	indm=i+ maxntheta*(*nphi-2)+ maxntheta* maxnphi*(*current_detector);
	rdphIP4dt=(int_psi4r_dt[indp]-int_psi4r_dt[indm])*idph;
	idphIP4dt=(int_psi4i_dt[indp]-int_psi4i_dt[indm])*idph;
      }
      dPrdp=rdphIP4dt;
      dPidp=idphIP4dt;

      // d_th I(P4)dt
      if (i>0 && i<*ntheta-1) {
	indp=(i+1)+ maxntheta*(j)+ maxntheta* maxnphi*(*current_detector);
	indm=(i-1)+ maxntheta*(j)+ maxntheta* maxnphi*(*current_detector);
	rdthIP4dt=(int_psi4r_dt[indp]-int_psi4r_dt[indm])*i2dth;
	idthIP4dt=(int_psi4i_dt[indp]-int_psi4i_dt[indm])*i2dth;
      }
      else if (i==0) {
	indp=1+ maxntheta*(j)+ maxntheta* maxnphi*(*current_detector); 
	indm=0+ maxntheta*(j)+ maxntheta* maxnphi*(*current_detector);
	rdthIP4dt=(int_psi4r_dt[indp]-int_psi4r_dt[indm])*idth;
	idthIP4dt=(int_psi4i_dt[indp]-int_psi4i_dt[indm])*idth;
      }
      else if (i==*ntheta-1) {
	indp=(*ntheta-1)+ maxntheta*(j)+ maxntheta* maxnphi*(*current_detector); 
	indm=(*ntheta-2)+ maxntheta*(j)+ maxntheta* maxnphi*(*current_detector);
	rdthIP4dt=(int_psi4r_dt[indp]-int_psi4r_dt[indm])*idth;
	idthIP4dt=(int_psi4i_dt[indp]-int_psi4i_dt[indm])*idth;
      }
      
      dPrdt=rdthIP4dt;
      dPidt=idthIP4dt;

      // complex conjugate of \int\int Psi_4 dtdt
      I2Pr= int2_psi4r_dt2[ind3d];
      I2Pi=-int2_psi4i_dt2[ind3d];

      f1=-sinp*dPrdt-cosp*cott*dPrdp-2*cosp*csct*iP4;
      f2=-sinp*dPidt-cosp*cott*dPidp+2*cosp*csct*rP4;
      dJxdt += -(f1*I2Pr-f2*I2Pi)*dOmega;

      f1=cosp*dPrdt-sinp*cott*dPrdp-2*sinp*csct*iP4;
      f2=cosp*dPidt-sinp*cott*dPidp+2*sinp*csct*rP4;
      dJydt += -(f1*I2Pr-f2*I2Pi)*dOmega;

      f1=dPrdp;
      f2=dPidp;
      dJzdt += -(f1*I2Pr-f2*I2Pi)*dOmega;
    }
  }
  // finish integrals
  const CCTK_REAL prefac=rad*rad/(4.*pi);

  dEdt =prefac*dEdt;
  dPxdt=prefac*dPxdt;
  dPydt=prefac*dPydt;
  dPzdt=prefac*dPzdt;
  dJxdt=prefac*dJxdt;
  dJydt=prefac*dJydt;
  dJzdt=prefac*dJzdt;

  const CCTK_INT det=*current_detector;

  // E(t)
  energy[det]=energy[det] + deltat * dEdt;

  // P^i(t)
  Px[det]=Px[det] + deltat * dPxdt;
  Py[det]=Py[det] + deltat * dPydt;
  Pz[det]=Pz[det] + deltat * dPzdt;

  // J^i
  Jx[det]=Jx[det] + deltat * dJxdt;
  Jy[det]=Jy[det] + deltat * dJydt;
  Jz[det]=Jz[det] + deltat * dJzdt;

  // IO
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"E=%g P^i=(%g,%g,%g) J^i=(%g,%g,%g)",
               energy[det],Px[det],Py[det],Pz[det],Jx[det],Jy[det],Jz[det]);
  }
  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"dE/dt=%g dP^i/dt=(%g,%g,%g) dJ^a/dt=(%g,%g,%g)",
               dEdt,dPxdt,dPydt,dPzdt,dJxdt,dJydt,dJzdt);
  }

  // update time
  lastt[*current_detector]=cctk_time;

  // write data only on CPU 0
  const CCTK_INT myproc=CCTK_MyProc(cctkGH);
  if (myproc==0) {
    CCTK_INT ierr=psi4analysis_write_data(CCTK_PASS_CTOC,
                                          det,
                                          rad,
                                          dEdt,
                                          dPxdt,dPydt,dPzdt,
					  dJxdt,dJydt,dJzdt);
    if (ierr<0) {
      CCTK_WARN(1,"data writing failed");
    }
  }
}



static int psi4analysis_write_data(CCTK_ARGUMENTS,
                                   CCTK_INT det,
                                   CCTK_REAL rad,
                                   CCTK_REAL dEdt,
                                   CCTK_REAL dPxdt, CCTK_REAL dPydt, CCTK_REAL dPzdt,
				   CCTK_REAL dJxdt, CCTK_REAL dJydt, CCTK_REAL dJzdt)
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
  if (*psi4analysis_files_exist>0) {
    sprintf(fmode,"a");
  }
  else {
    sprintf(fmode,"w");
  }

  // filename
  if (name_detector_by_radius) {
      sprintf(varname, "psi4analysis_r%.2f",rad);
  }
  else {
      sprintf(varname,"psi4analysis_det_%d",det);
  }

  filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                              strlen (file_extension) +2);
  sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_ylm_mode: Could not open scalar output file '%s'",
                filename);
    return -1;
  }

  // write header on startup
  if (*psi4analysis_files_exist<=0) {
    fprintf(file,"# Psi4 Analysis Quantities\n");
    fprintf(file,"# det no.=%d rdet=%g\n",det,rad);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:t 2:E 3:Px 4:Py 5:Pz 6:dE/dt 7:dPx/dt 8:dPy/dt 9:dPz/dt 10:dJ^x/dt 11:dJ^y/dt 12:dJ^z/dt 13:J^x 14:J^y 15:J^z\n");
  }

  // write data
  sprintf (format_str_real,
           "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
           out_format,out_format,out_format,out_format,out_format,out_format,
           out_format,out_format,out_format,out_format,out_format,out_format,
	   out_format,out_format,out_format);
  fprintf(file, format_str_real, cctk_time, energy[det],
          Px[det], Py[det], Pz[det],
          dEdt,
	  dPxdt,dPydt,dPzdt,
	  dJxdt,dJydt,dJzdt,
	  Jx[det],Jy[det],Jz[det]);

  fclose(file); 
  free(filename); 
  if (det==number_of_detectors-1) { // we have written all detectors
    *psi4analysis_files_exist=1;
  }
  return 1;
}
