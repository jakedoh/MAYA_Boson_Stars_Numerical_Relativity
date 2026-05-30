#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 


static int scalaranalysis_write_data(CCTK_ARGUMENTS,
                                   CCTK_INT det,
                                   CCTK_REAL rad,
                                   CCTK_REAL dEdt,
                                   CCTK_REAL dPxdt, CCTK_REAL dPydt, CCTK_REAL dPzdt, CCTK_REAL dJzdt);


void scalarAnalysis_Compute_EPRad (CCTK_ARGUMENTS)
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
  const CCTK_REAL dtp=dth*dph;
  
  CCTK_REAL th=0,ph=0,eTtr=0,eTrphi=0,eTrx=0,eTry=0,eTrz=0;
  CCTK_REAL sint=0,cost=0,sinp=0,cosp=0,cott=0,csct=0;
  CCTK_REAL nx=0,ny=0,nz=0;
  CCTK_REAL nphix=0,nphiy=0,nphiz=0;
  CCTK_REAL iwtheta=0,iwphi=0,intweight=0;

  // surface integrals
  CCTK_REAL dEdt=0;
  CCTK_REAL dPxdt=0,dPydt=0,dPzdt=0;
  CCTK_REAL dJzdt=0;
  CCTK_INT ind3d;
  CCTK_INT ind2d=0;

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
      ind2d=i+ maxntheta*j;
      ind3d=i+ maxntheta*j+ maxntheta* maxnphi*(*current_detector);

      nx=cosp*sint;
      ny=sinp*sint;
      nz=     cost;

      eTtr=nx*Tmunu_2d_eTtx[ind2d]+ny*Tmunu_2d_eTty[ind2d]+nz*Tmunu_2d_eTtz[ind2d];
      // linear momentum
      // radiated quantities
      // energy
      dEdt+=-eTtr*dOmega;

      dPxdt+=(nx*Tmunu_2d_eTxx[ind2d]+ny*Tmunu_2d_eTxy[ind2d]+nz*Tmunu_2d_eTxz[ind2d])*dOmega;
      dPydt+=(nx*Tmunu_2d_eTxy[ind2d]+ny*Tmunu_2d_eTyy[ind2d]+nz*Tmunu_2d_eTyz[ind2d])*dOmega;
      dPzdt+=(nx*Tmunu_2d_eTxz[ind2d]+ny*Tmunu_2d_eTyz[ind2d]+nz*Tmunu_2d_eTzz[ind2d])*dOmega;

      nphix= -rad*sint*sinp;
      nphiy=  rad*sint*cosp;
      nphiz=              0;

      eTrx=nx*Tmunu_2d_eTxx[ind2d]+ny*Tmunu_2d_eTxy[ind2d]+nz*Tmunu_2d_eTxz[ind2d];
      eTry=nx*Tmunu_2d_eTxy[ind2d]+ny*Tmunu_2d_eTyy[ind2d]+nz*Tmunu_2d_eTyz[ind2d];
      eTrz=nx*Tmunu_2d_eTxz[ind2d]+ny*Tmunu_2d_eTyz[ind2d]+nz*Tmunu_2d_eTzz[ind2d];

      eTrphi=nphix*eTrx+nphiy*eTry+nphiz*eTrz;
      dJzdt+=eTrphi*dOmega;
    }
  }
  // finish integrals
  const CCTK_REAL prefac=rad*rad;

  dEdt =prefac*dEdt;
  dPxdt=prefac*dPxdt;
  dPydt=prefac*dPydt;
  dPzdt=prefac*dPzdt;
  dJzdt=prefac*dJzdt;

  const CCTK_INT det=*current_detector;

  // E(t)
  energy[det]=energy[det] + deltat * dEdt;

  // P^i(t)
  Px[det]=Px[det] + deltat * dPxdt;
  Py[det]=Py[det] + deltat * dPydt;
  Pz[det]=Pz[det] + deltat * dPzdt;

  Jz[det]=Jz[det] + deltat * dJzdt;
  // IO
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"E=%g P^i=(%g,%g,%g) Jz=%g",
               energy[det],Px[det],Py[det],Pz[det],Jz[det]);
  }
  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"dE/dt=%g dP^i/dt=(%g,%g,%g)",
               dEdt,dPxdt,dPydt,dPzdt);
  }

  // update time
  lastt[*current_detector]=cctk_time;

  // write data only on CPU 0
  const CCTK_INT myproc=CCTK_MyProc(cctkGH);
  if (myproc==0) {
    CCTK_INT ierr=scalaranalysis_write_data(CCTK_PASS_CTOC,
                                          det,
                                          rad,
                                          dEdt,
                                          dPxdt,dPydt,dPzdt,dJzdt);
    if (ierr<0) {
      CCTK_WARN(1,"data writing failed");
    }
  }
}



static int scalaranalysis_write_data(CCTK_ARGUMENTS,
                                   CCTK_INT det,
                                   CCTK_REAL rad,
                                   CCTK_REAL dEdt,
                                   CCTK_REAL dPxdt, CCTK_REAL dPydt, CCTK_REAL dPzdt, CCTK_REAL dJzdt
				  )
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char fmode[1];
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  // file mode: append if already written
  if (*scalaranalysis_files_exist>0) {
    sprintf(fmode,"a");
  }
  else {
    sprintf(fmode,"w");
  }

  // filename
  if (name_detector_by_radius) {
      sprintf(varname, "scalaranalysis_r%.2f",rad);
  }
  else {
      sprintf(varname,"scalaranalysis_det_%d",det);
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
  if (*scalaranalysis_files_exist<=0) {
    fprintf(file,"# scalar Analysis Quantities\n");
    fprintf(file,"# det no.=%d rdet=%g\n",det,rad);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:t 2:E 3:Px 4:Py 5:Pz 6:Jz 7:dE/dt 8:dPx/dt 9:dPy/dt 10:dPz/dt 11:dJz/dt\n");
  }

  // write data
  sprintf (format_str_real,
           "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
           out_format,out_format,out_format,out_format,out_format,out_format,
           out_format,out_format,out_format,out_format,out_format);
  fprintf(file, format_str_real, cctk_time, energy[det],
          Px[det], Py[det], Pz[det], Jz[det],
          dEdt, dPxdt, dPydt, dPzdt, dJzdt);

  fclose(file); 
  free(filename); 
  if (det==number_of_detectors-1) { // we have written all detectors
    *scalaranalysis_files_exist=1;
  }
  return 1;
}
