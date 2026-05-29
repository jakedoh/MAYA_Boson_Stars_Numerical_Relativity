#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"


int write_moncriefq_mode (const CCTK_REAL time,
                          const CCTK_INT wavextr_files_exist,
                          const CCTK_INT first_detector_this_time,
                          const CCTK_INT idet,
                          const CCTK_REAL rad,
                          const CCTK_INT il,
                          const CCTK_INT im,
                          const CCTK_REAL qevenr, // real part
                          const CCTK_REAL qeveni, // imag part
                          const CCTK_REAL qoddr,  // real part
                          const CCTK_REAL qoddi,  // imaginary part
                          const CCTK_REAL rschw,
                          const CCTK_REAL mschw)
{
  DECLARE_CCTK_PARAMETERS;

  char fmode[1];
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[1024]; // XXX fixed size
  FILE *file;
  CCTK_REAL phase,amp;


  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "in write_moncriefq_mode");
  }

  // file mode: append if already written
  if (wavextr_files_exist>0 ||
      (first_detector_this_time<=0 && write_single_file_per_mode>0)) {
    sprintf(fmode,"a");
  }
  else {
    sprintf(fmode,"w");
  }

  // filename
  if (write_single_file_per_mode) {
    sprintf(varname, "MoncriefQ_l%d_m%d",il,im);
  }
  else {
    sprintf(varname, "MoncriefQ_l%d_m%d_r%4.2f",il,im,rad);
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

  // file header
  if (wavextr_files_exist<=0) {
    if (write_single_file_per_mode) {
      fprintf(file,"# MoncriefQ\n");
      fprintf(file,"#    l=%d m=%d\n",il,im);
      fprintf(file,"# 1:r 2:re(Qeven) 3:im(Qeven) 4:re(Qodd) 5:im(Qodd) 6:r_schw 7:m_schw\n");
    }
    else {
      fprintf(file,"# MoncriefQ\n");
      fprintf(file,"#    l=%d m=%d r=%f ndet=%d\n",il,im,rad,idet);
      fprintf(file,"# 1:t 2:re(Qeven) 3:im(Qeven) 4:re(Qodd) 5:im(Qodd) 6:r_schw 7:m_schw\n");
    }
  }

  // single file per detector or all detectors into one file
  sprintf (format_str_real, "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n",
           out_format,out_format,out_format,out_format,out_format,
           out_format,out_format);

  if (write_single_file_per_mode) {
    if (first_detector_this_time>0) {
      fprintf(file, "\n\n#Time = %g\n", time);
    }
    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "write rad=%g rQ+=%g iQ+=%g rQx=%g iQx=%g r_s=%g, m_s=%g into %s",
                 rad,qevenr,qeveni,qoddr,qoddi,rschw,mschw,filename);
    }
    fprintf(file, format_str_real, rad, qevenr, qeveni, qoddr, qoddi, rschw, mschw);
  }
  else {
    fprintf(file, format_str_real, time, qevenr, qeveni,
                                   qoddr, qoddi, rschw, mschw);
  }

  fclose(file);
  free(filename);

  return 1;
}



void CCTK_FCALL CCTK_FNAME(write_moncriefq_mode)
  (CCTK_INT *ierr,
   const CCTK_REAL *time,
   const CCTK_INT *wavextr_files_exist,
   const CCTK_INT *first_detector_this_time,
   const CCTK_INT *idet,
   const CCTK_REAL *rad,
   const CCTK_INT *il,
   const CCTK_INT *im,
   const CCTK_REAL *qevenr, // real part
   const CCTK_REAL *qeveni, // imag part
   const CCTK_REAL *qoddr,  // real part
   const CCTK_REAL *qoddi,  // imaginary part
   const CCTK_REAL *rschw,
   const CCTK_REAL *mschw)
{
  DECLARE_CCTK_PARAMETERS;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING,"in write_moncriefq_mode fortran wrapper at T=%f",
               *time);
    CCTK_VInfo(CCTK_THORNSTRING,"  time=%f, wavestr_files_exist=%d, first_detector_this_time=%d",
               *time,*wavextr_files_exist,*first_detector_this_time);
    CCTK_VInfo(CCTK_THORNSTRING,"  idet=%d,rad=%f,il=%d,im=%d",
               *idet,*rad,*il,*im);
    CCTK_VInfo(CCTK_THORNSTRING,"  rQ+=%g, iQ+=%g, rQx=%g iQx=%g",
               *qevenr, *qeveni, *qoddr, *qoddi);
  }

  *ierr=write_moncriefq_mode(*time,
                       *wavextr_files_exist,
                       *first_detector_this_time,
                       *idet,
                       *rad,
                       *il,
                       *im,
                       *qevenr, // real part
                       *qeveni, // imag part
                       *qoddr,  // real part
                       *qoddi,  // imaginary part
                       *rschw,
                       *mschw);
}



int write_radiated_energy_momentum(CCTK_INT cdt, CCTK_INT create_file,
                                   CCTK_REAL rad, CCTK_REAL time,
                                   CCTK_REAL E, CCTK_REAL dEdt,
                                   CCTK_REAL Px, CCTK_REAL Py, CCTK_REAL Pz,
                                   CCTK_REAL dPxdt, CCTK_REAL dPydt,
                                   CCTK_REAL dPzdt)
{
  DECLARE_CCTK_PARAMETERS;

  char fmode[1];
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  FILE *file;

  if (verbose>3) {
    CCTK_INFO("in write_radiated_energy_momentum (c function)");
    CCTK_VInfo(CCTK_THORNSTRING,"det=%d create_file=%d rad=%g E=%g dEdt=%g",
               cdt,create_file,rad,E,dEdt);
    CCTK_VInfo(CCTK_THORNSTRING,"P^i=[%g,%g,%g] dP^i/dt=[%g,%g,%g]",
               Px,Py,Pz,dPxdt,dPydt,dPzdt);
  }

  // file mode: append if already written
  if (create_file<=0) {
    sprintf(fmode,"a");
  }
  else {
    sprintf(fmode,"w");
  }

  // filename
  sprintf(varname, "Zerilli_radEP_det_%d",cdt);

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
  if (create_file>0) {
    fprintf(file,"# Zerilli radiated energy and linear momentum Quantities\n");
    fprintf(file,"# det. no. %d rad=%g\n",cdt,rad);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:t 2:E 3:Px 4:Py 5:Pz 6:dE/dt 7:dPx/dt 8:dPy/dt 9:dPz/dt\n");
  }

  // write data
  sprintf (format_str_real,
           "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
           out_format, out_format,out_format,out_format,out_format,out_format,
           out_format,out_format,out_format);
  fprintf(file, format_str_real, time, E,
          Px, Py, Pz,
          dEdt,dPxdt,dPydt,dPzdt);

  fclose(file); 
  free(filename); 

  return 1;
}

void CCTK_FCALL CCTK_FNAME(write_radiated_energy_momentum)
  (CCTK_INT *ierr,
   CCTK_INT *cdt,
   CCTK_INT *create_file,
   CCTK_REAL *rad,
   CCTK_REAL *time,
   CCTK_REAL *E,
   CCTK_REAL *dEdt,
   CCTK_REAL *Px,
   CCTK_REAL *Py,
   CCTK_REAL *Pz,
   CCTK_REAL *dPxdt,
   CCTK_REAL *dPydt,
   CCTK_REAL *dPzdt)
{
  DECLARE_CCTK_PARAMETERS;

  if (verbose>3) {
    CCTK_INFO("in write_radiated_energy_momentum (fortran wrapper)");
    CCTK_VInfo(CCTK_THORNSTRING,"det=%d creaf=%d E=%g dEdt=%g",
               *cdt,*create_file,*E,*dEdt);
    CCTK_VInfo(CCTK_THORNSTRING,"P^i=[%g,%g,%g] dP^i/dt=[%g,%g,%g]",
               *Px,*Py,*Pz,*dPxdt,*dPydt,*dPzdt);
  }
   
  *ierr = write_radiated_energy_momentum(*cdt,*create_file,*rad,*time,
                                         *E,*dEdt,*Px,*Py,*Pz,
                                         *dPxdt,*dPydt,*dPzdt);

}
