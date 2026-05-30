 /*@@
   @file      IO.c
   @date      Mon Nov 25 13:19:46 2002
   @author    Frank Herrmann
   @desc 
              Write output data to file.
   @enddesc 
   @version  $Id: IO.c,v 1.14 2004/10/06 16:08:04 herrmann Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

static const char *rcsid = "$Header: /numrelcvs/HerrmannCVS/Ylm_Decomp/src/IO.c,v 1.14 2004/10/06 16:08:04 herrmann Exp $";

CCTK_FILEVERSION(HerrmannCVS_Ylm_Decomp_Write_Data_c);

int write_ylm_mode (const CCTK_REAL time,
                    const CCTK_INT ylm_files_exist,
                    const CCTK_INT cmplx_var, // complex variable?
                    const CCTK_INT first_detector_this_time,
                    const CCTK_INT idet,
                    const CCTK_REAL rad,
                    const CCTK_INT vind,
                    const CCTK_INT il,
                    const CCTK_INT im,
                    const CCTK_REAL rval, // real part
                    const CCTK_REAL ival)  // imaginary part
{
  DECLARE_CCTK_PARAMETERS;

  char *filename = NULL;
  char varname[1024]; // XXX fixed size
  const char *file_extension = ".asc";
  char format_str_real[1024]; // XXX fixed size
  FILE *file = NULL;

  char *var_name = CCTK_FullName(vind); // Remember to free this

  //----------------------------------------------


  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "in write_ylm_mode");
  }


  // filename
  if (write_single_file_per_mode) {
    sprintf(varname, "Ylm_%s_l%d_m%d",
            var_name,il,im);
  }
  else {
    if (name_detector_by_radius) {
      sprintf(varname, "Ylm_%s_l%d_m%d_r%4.2f",
              var_name,il,im,rad);
    }
    else {
      sprintf(varname, "Ylm_%s_l%d_m%d_n%04d",
              var_name,il,im,idet);
    }
  }
  filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                              strlen (file_extension) +2);
  sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

  // open file
  file = fopen (filename, 
                (ylm_files_exist > 0 ||
                 (first_detector_this_time<=0 && write_single_file_per_mode>0)
                 ? "a" : "w" ));  // file mode: append if already written

  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_ylm_mode: Could not open scalar output file '%s'",
                filename);
    return -1;
  }

  // file header
  if (ylm_files_exist<=0) {
    if (write_single_file_per_mode) {
      fprintf(file,"# Ylm_Decomp\n");
      fprintf(file,"#    %s l=%d m=%d\n",var_name,il,im);
      if (cmplx_var>0) {
        fprintf(file,"# 1:t 2:re(y) 3:im(y)\n");
      }
    }
    else {
      fprintf(file,"# Ylm_Decomp\n");
      fprintf(file,"#    %s l=%d m=%d r=%f\n",
              var_name,il,im,rad);
      if (cmplx_var>0) {
        fprintf(file,"# 1:t 2:re(y) 3:im(y)\n");
      }
    }
  }

  // single file per detector or all detectors into one file
  if (cmplx_var>0) {
    sprintf (format_str_real, "%%%s\t%%%s\t%%%s\n",
             out_format,out_format,out_format);
  }
  else {
    sprintf (format_str_real, "%%%s\t%%%s\n", out_format, out_format);
  }
  if (write_single_file_per_mode) {
    if (first_detector_this_time>0) {
      fprintf(file, "\n\n#Time = %g\n", time);
    }
    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,"write rad=%g rval=%g ival=%g into %s",
                 rad,rval,ival,filename);
    }
    if (cmplx_var>0) {
      fprintf(file, format_str_real, rad, rval, ival);
    }
    else {
      fprintf(file, format_str_real, rad, rval);
    }
  }
  else {
    if (cmplx_var>0) {
      fprintf(file, format_str_real, time, rval, ival);
    }
    else {
      fprintf(file, format_str_real, time, rval);
    }
  }

  fclose(file);
  free(filename);
  free(var_name);

  return 1;
}

void CCTK_FCALL CCTK_FNAME(write_ylm_mode)
  (CCTK_INT *ierr,
   const CCTK_REAL *time,
   const CCTK_INT *ylm_files_exist,
   const CCTK_INT *cmplx_var,
   const CCTK_INT *first_detector_this_time,
   const CCTK_INT *idet,
   const CCTK_REAL *rad,
   const CCTK_INT *vind,
   const CCTK_INT *il,
   const CCTK_INT *im,
   const CCTK_REAL *rval,
   const CCTK_REAL *ival)
{
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING,"in write_ylm_mode fortran wrapper at T=%f",
               *time);
    CCTK_VInfo(CCTK_THORNSTRING,"  time=%f, ylm_files_exist=%d, first_detector_this_time=%d, complex_var=%d",
               *time,*ylm_files_exist,*first_detector_this_time,*cmplx_var);
    CCTK_VInfo(CCTK_THORNSTRING,"  idet=%d,rad=%f,vind=%d,il=%d,im=%d,rval=%g,ival=%g",
               *idet,*rad,*vind,*il,*im,*rval,*ival);
  }

  *ierr=write_ylm_mode(*time, *ylm_files_exist, *cmplx_var,
                       *first_detector_this_time,
                       *idet,*rad,*vind,*il,*im,*rval,*ival);  
}








// write out 2D surface quantities if required
// integrand characterises data: 0 means f*Ylm, 1 means f
//   2: spinw -2 real
//   3: spinw -2 imag
//   4: spinw 0 real
//   5: spinw 0 imag
//   6: spinw 2 real
//   7: spinw 2 imag
int write_surface_quantities_mode
  (const CCTK_REAL time,
   const CCTK_INT ylm_files_exist,
   const CCTK_INT first_detector_this_time,
   const CCTK_INT idet,
   const CCTK_REAL rad,
   const CCTK_INT vind,
   const CCTK_INT il,
   const CCTK_INT im,
   const CCTK_REAL *vals,
   const CCTK_INT integrand,
   const CCTK_INT lntheta,
   const CCTK_INT lnphi)
{
  DECLARE_CCTK_PARAMETERS;

  char fmode[1];
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[1024]; // XXX fixed size
  FILE *file;
  char *var_name = CCTK_FullName(vind); // Remember to free this

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "in write_ylm_mode");
  }

  // file mode: append if already written
  if (ylm_files_exist>0 || first_detector_this_time<=0) {
    sprintf(fmode,"a");
  }
  else {
    sprintf(fmode,"w");
  }

  // filename depends on what we write out
  switch (integrand)
    {
    case 1 :
      sprintf(varname, "Ylm_%s_l%d_m%d_r%4.2f_2d_integrand",
              var_name,il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing 2D integrand");
      }
      break;

    case 0 :
      sprintf(varname, "Ylm_%s_l%d_m%d_r%4.2f_2d_variable",
              var_name,il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing 2D variable");
      }
      break;

    case 2 :
      sprintf(varname, "rs-2Ylm_l%d_m%d_r%4.2f",il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing spinw -2 Ylm real");
      }
      break;

    case 3 :
      sprintf(varname, "is-2Ylm_l%d_m%d_r%4.2f",il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing spinw -2 Ylm imag");
      }
      break;

    case 4 :
      sprintf(varname, "rYlm_l%d_m%d_r%4.2f",il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing spinw 0 Ylm real");
      }
      break;

    case 5 :
      sprintf(varname, "iYlm_l%d_m%d_r%4.2f",il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing spinw 0 Ylm imag");
      }
      break;

    case 6 :
      sprintf(varname, "rs2Ylm_l%d_m%d_r%4.2f",il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing spinw 2 Ylm real");
      }
      break;

    case 7 :
      sprintf(varname, "is2Ylm_l%d_m%d_r%4.2f",il,im,rad);
      if (verbose>3) {
        CCTK_INFO("writing spinw 2 Ylm imag");
      }
      break;

    default:
      if (verbose>2) {
        CCTK_WARN(1,"unclear mode for write_surface_q...");
      }
      return -1;
      break;
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
  if (ylm_files_exist<=0) {
    fprintf(file,"# Ylm_Decomp 2D surface\n");
    fprintf(file,"#    %s l=%d m=%d\n",var_name,il,im);
    fprintf(file,"#    theta phi value\n");
  }

  sprintf (format_str_real, "%%%s\t%%%s\t%%%s\n", out_format, out_format, out_format);
  if (first_detector_this_time>0) {
    fprintf(file, "\n\n#Time = %g\n", time);
  }

  // dump it all out : theta in [0,pi], phi in [0,2 pi]
  const CCTK_REAL pi=acos(-1.0);
  const CCTK_REAL dth=pi/(lntheta-1);
  const CCTK_REAL dph=2.*pi/(lnphi-1);
  CCTK_REAL th,ph;
  CCTK_INT ind;
  for (int j=0;j<lnphi;j++) {
    ph=(j+0.5)*dph;
    for (int i=0;i<lntheta;i++) {
      th=(i+0.5)*dth; // stagger the axis
      ind=i+lntheta*j; // fortran index confusion
      fprintf(file, format_str_real, th,ph,vals[ind]);
    }
    if (j!=lnphi-1) {
      fprintf(file,"\n"); // new line indicator for gnuplot
    }
  }

  fclose(file);
  free(filename);
  free(var_name);

  return 1;
}

void CCTK_FCALL CCTK_FNAME(write_surface_quantities_mode)
  (CCTK_INT *ierr,
   const CCTK_REAL *time,
   const CCTK_INT *ylm_files_exist,
   const CCTK_INT *first_detector_this_time,
   const CCTK_INT *idet,
   const CCTK_REAL *rad,
   const CCTK_INT *vind,
   const CCTK_INT *il,
   const CCTK_INT *im,
   const CCTK_REAL *vals,
   const CCTK_INT *integrand,
   const CCTK_INT *lntheta,
   const CCTK_INT *lnphi)
{
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "in write_surface_quantities_mode fortran wrapper at T=%f",
               *time);
    CCTK_VInfo(CCTK_THORNSTRING,
               "  time=%f, ylm_files_exist=%d, first_detector_this_time=%d, integrand=%d",
               *time,*ylm_files_exist,*first_detector_this_time,*integrand);
    CCTK_VInfo(CCTK_THORNSTRING,
               "  idet=%d,rad=%f,vind=%d,il=%d,im=%d,lntheta=%d,lnphi=%d,val0=%g",
               *idet,*rad,*vind,*il,*im,*lntheta,*lnphi,vals[0]);
  }

  *ierr=write_surface_quantities_mode(*time, *ylm_files_exist,
                                      *first_detector_this_time,
                                      *idet,*rad,*vind,*il,*im,
                                      vals,*integrand,*lntheta,*lnphi);  
}
