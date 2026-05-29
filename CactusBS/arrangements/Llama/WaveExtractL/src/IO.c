
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

 /*@@
   @file      IO.c
   @date      Mon Nov 25 13:19:46 2002
   @author    Frank Herrmann
   @desc 
              Write output data to file.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "cctk_Timers.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extractGH.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"


void WavExtrL_WriteData(CCTK_ARGUMENTS);

void WavExtrL_WriteLotsOfData(CCTK_ARGUMENTS);

int WavExtrL_TimeForOutput (const cGH *GH, int vindex);

static void WavExtrL_WriteTimer(cTimerData *info);

static void *WavExtrL_SetupGH (tFleshConfig *config, int conv_level, cGH *GH);

static int WavExtrL_WriteScalar(const cGH *GH,
                        CCTK_INT vindex,
                        CCTK_REAL  *value,
                        CCTK_INT   max,
                        const char *alias,
                        const CCTK_REAL  *detector_radii);

static FILE *WavExtrL_OpenFile (const cGH *GH,
                               int vindex,
                               const char *filename,
                               const char *slicename,
                               const char *description,
                               CCTK_INT   number_sets,
                               const char *aliasname,
                               const CCTK_REAL  *detector_radii);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

 /*@@
   @routine   WavExtrL_Startup
   @date      Mon 16th December 2002
   @author    Gabrielle Allen
   @desc
              The startup registration routine for Extract.
              Registers WaveExtract as an IO Method and provide the
              only necessary method TimeForOutput
   @enddesc
   @calls     CCTK_RegisterGHExtensionSetupGH
@@*/
int WavExtrL_Startup (void)
{
  int ierr;

  if (CCTK_GHExtensionHandle ("IO") < 0)
  {
    CCTK_WARN (1, "Thorn IOUtil was not activated. "
                  "No WaveExtractL I/O methods will be enabled.");
    return 0;
  }

  CCTK_RegisterGHExtensionSetupGH(CCTK_RegisterGHExtension("WaveExtractL"),
                                   WavExtrL_SetupGH);

  ierr = CCTK_TimerCreate("MoncriefQ");
  if (ierr < 0)
  {
    CCTK_WARN(1,"could not create timer MoncriefQ");
  }

  ierr =  CCTK_TimerCreate("Interpolation");
  if (ierr < 0)
  {
    CCTK_WARN(1,"could not create timer Interpolation");
  }
   
  ierr =  CCTK_TimerCreate("Schwarzschild");
  if (ierr < 0)
  {
    CCTK_WARN(1,"could not create timer Schwarzschild");
  }

  return 0;
}


 /*@@
   @routine    WavExtrL_WriteData
   @date       Mon Nov 25 13:20:06 2002
   @author     Frank Herrmann
   @desc 
               Write output data to file.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
void WavExtrL_WriteData(CCTK_ARGUMENTS)
{
  char string_var[256],radiusname[256];
  CCTK_INT vindex_rsch;
  CCTK_INT vindex_mass;
  CCTK_INT vindex_Qeven_Re;
  CCTK_INT vindex_Qodd_Re;
  CCTK_INT vindex_Qeven_Im;
  CCTK_INT vindex_Qodd_Im;
  CCTK_INT il,im;
  CCTK_INT index;
  CCTK_REAL *data_Qeven_Re;
  CCTK_REAL *data_Qodd_Re;
  CCTK_REAL *data_Qeven_Im;
  CCTK_REAL *data_Qodd_Im;
  CCTK_REAL val;

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  if (verbose>2)
    CCTK_INFO("Writing output data to individual files per detector");

  if (*do_nothing==1)
    return;

  if (cctk_iteration != 0) {
    if ((cctk_iteration % my_out_every_det[*current_detector])!=0) {
      if (verbose>2) CCTK_INFO("No time for this detector");
      return;
    }
  }

  if (calc_when_necessary == 1)
  {
    if (cctk_time < *current_detector_radius-50)
    {
      if (verbose>2) CCTK_INFO("No time for this detector");
      return;
    } 
    if (CCTK_IsFunctionAliased("MergerHandler_WeHaveMerger"))
    {
      if (MergerHandler_WeHaveMerger() == 1)
      {
        if (cctk_time > MergerHandler_MergerTime()+*current_detector_radius+ringdown_margin)
	{
          if (verbose>2) CCTK_INFO("No time for this detector");
          return;
        }
      }
    }
  } 


  vindex_rsch  = CCTK_VarIndex("WaveExtractL::Schwarzschild_Radius");
  vindex_mass  = CCTK_VarIndex("WaveExtractL::Schwarzschild_Mass");
  vindex_Qeven_Re = CCTK_VarIndex("WaveExtractL::Qeven_Re");
  vindex_Qodd_Re  = CCTK_VarIndex("WaveExtractL::Qodd_Re");
  vindex_Qeven_Im = CCTK_VarIndex("WaveExtractL::Qeven_Im");
  vindex_Qodd_Im  = CCTK_VarIndex("WaveExtractL::Qodd_Im");

  if(surface_index[*current_detector]>=0) {
    snprintf(radiusname, sizeof radiusname,
             "Detector_Number_%d",(int)*current_detector);
  }
  else {
    snprintf(radiusname, sizeof radiusname,
             "Detector_Radius_%4.2f",(double)*current_detector_radius);
  }

  if(CCTK_MyProc (cctkGH) == 0)
  {
    /* Schw Mass indicator */
    snprintf(string_var, sizeof string_var,
             "Schwarzschild_Mass_%s",radiusname);
    WavExtrL_WriteScalar(cctkGH,vindex_mass,Schwarzschild_Mass,
                        -1,string_var,NULL);

    /* Schw Radius */
    snprintf(string_var, sizeof string_var,
             "Schwarzschild_Radius_%s",radiusname);
    WavExtrL_WriteScalar(cctkGH,vindex_rsch,Schwarzschild_Radius,
                        -1,string_var,NULL);

    /* write information to screen about the results */
    if(verbose>1)
    {
      CCTK_VInfo(CCTK_THORNSTRING,"Detector No. %d",(int)*current_detector);
      CCTK_VInfo(CCTK_THORNSTRING,
              "  Coordinate Radius    = %4.3f",(double)*current_detector_radius);
      CCTK_VInfo(CCTK_THORNSTRING,
              "  Schwarzschild Radius = %2.3f",(double)*Schwarzschild_Radius);
      CCTK_VInfo(CCTK_THORNSTRING,
              "  Schwarzschild Mass   = %2.3f",(double)*Schwarzschild_Mass);
    }

    data_Qeven_Re=(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qeven_Re);
    data_Qodd_Re=(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qodd_Re);
    data_Qeven_Im=(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qeven_Im);
    data_Qodd_Im=(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qodd_Im);

    /* loop over l,m modes */
    for (il=*l_min; il<=*l_max; il=il+*l_step) {
      for (im=*m_min; im<=*m_max; im=im+*m_step) {
        if ((im<-il) || (im>+il)) continue;
        /* Qeven */
        index=(il-1)+im*(*l_max);
        snprintf(string_var, sizeof string_var,
                 "Qeven_Re_%s_l%d_m%d",radiusname,(int)il,(int)im);
        val=data_Qeven_Re[index];
        WavExtrL_WriteScalar(cctkGH,vindex_Qeven_Re,
                         &val,-1,string_var,NULL);
        snprintf(string_var, sizeof string_var,
                 "Qeven_Im_%s_l%d_m%d",radiusname,(int)il,(int)im);
        val=data_Qeven_Im[index];
        WavExtrL_WriteScalar(cctkGH,vindex_Qeven_Im,
                         &val,-1,string_var,NULL);
        /* Qodd */
        snprintf(string_var, sizeof string_var,
                 "Qodd_Re_%s_l%d_m%d",radiusname,(int)il,(int)im);
        val=data_Qodd_Re[index];
        WavExtrL_WriteScalar(cctkGH,vindex_Qodd_Re,
                         &val,-1,string_var,NULL);
        snprintf(string_var, sizeof string_var,
                 "Qodd_Im_%s_l%d_m%d",radiusname,(int)il,(int)im);
        val=data_Qodd_Im[index];
        WavExtrL_WriteScalar(cctkGH,vindex_Qodd_Im,
                         &val,-1,string_var,NULL);
        if (verbose>3) {
          CCTK_VInfo(CCTK_THORNSTRING,
             "    IO index: %d (l,m)=(%d,%d) Re: Qeven: %16.8g Qodd: %16.8g",
                  (int)index,(int)il,(int)im,(double)data_Qeven_Re[index],(double)data_Qodd_Re[index]);
           CCTK_VInfo(CCTK_THORNSTRING,
             "    IO index: %d (l,m)=(%d,%d) Im: Qeven: %16.8g Qodd: %16.8g",
                  (int)index,(int)il,(int)im,(double)data_Qeven_Im[index],(double)data_Qodd_Im[index]);
        }
      }
    }
  }
}

/* for a lot of detectors we have a different output format */
void WavExtrL_WriteLotsOfData(CCTK_ARGUMENTS)
{
  char string_var[256];
  CCTK_INT vindex_rsch;
  CCTK_INT vindex_mass;
  CCTK_INT vindex_Qeven_Re;
  CCTK_INT vindex_Qodd_Re;
  CCTK_INT vindex_Qeven_Im;
  CCTK_INT vindex_Qodd_Im;
  CCTK_INT index;
  CCTK_INT il,im;
  CCTK_REAL *data_Qeven_Re;
  CCTK_REAL *data_Qodd_Re;
  CCTK_REAL *data_Qeven_Im;
  CCTK_REAL *data_Qodd_Im;
  CCTK_REAL *val_even_Re, *val_odd_Re;
  CCTK_REAL *val_even_Im, *val_odd_Im;
  CCTK_INT idet;

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  if (verbose>4)
    CCTK_INFO("Writing output data to many files");

  if (*do_nothing==1)
    return;
  vindex_rsch  = CCTK_VarIndex("WaveExtractL::Schw_Radii");
  vindex_mass  = CCTK_VarIndex("WaveExtractL::Schw_Masses");
  vindex_Qeven_Re = CCTK_VarIndex("WaveExtractL::Qeven_Re_Array");
  vindex_Qodd_Re  = CCTK_VarIndex("WaveExtractL::Qodd_Re_Array");
  vindex_Qeven_Im = CCTK_VarIndex("WaveExtractL::Qeven_Im_Array");
  vindex_Qodd_Im  = CCTK_VarIndex("WaveExtractL::Qodd_Im_Array");

  if (maximum_detector_number==1)
    return;

  if(CCTK_MyProc (cctkGH) == 0)
  {
    /* Schw Mass indicator */
    snprintf(string_var, sizeof string_var,
             "Schwarzschild_Masses");
    WavExtrL_WriteScalar(cctkGH,vindex_mass,Schw_Masses,
                        maximum_detector_number,
                        string_var,detector_radius);
    /* Schw Radius */
    snprintf(string_var, sizeof string_var,
             "Schwarzschild_Radii");
    WavExtrL_WriteScalar(cctkGH,vindex_rsch,Schw_Radii,
                        maximum_detector_number,
                        string_var,detector_radius);

    /* write information to screen about the results */
    if(verbose>1)
    {
      CCTK_VInfo(CCTK_THORNSTRING,
         "Writing data from %d detectors to disk",(int)maximum_detector_number);
      CCTK_VInfo(CCTK_THORNSTRING,"Last Detector data (No. %d)",
                 (int)*current_detector);
      CCTK_VInfo(CCTK_THORNSTRING,
              "  Coordinate Radius    = %4.3f",(double)*current_detector_radius);
      CCTK_VInfo(CCTK_THORNSTRING,
              "  Schwarzschild Radius = %2.3f",(double)*Schwarzschild_Radius);
      CCTK_VInfo(CCTK_THORNSTRING,
              "  Schwarzschild Mass   = %2.3f",(double)*Schwarzschild_Mass);
    }

    /* output of the Moncrief Qeven/Qodd */
    data_Qeven_Re=(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qeven_Re);
    data_Qodd_Re =(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qodd_Re);
    data_Qeven_Im=(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qeven_Im);
    data_Qodd_Im =(CCTK_REAL *) CCTK_VarDataPtrI(cctkGH,0,vindex_Qodd_Im);
    val_even_Re=(CCTK_REAL *) malloc(*max_det_no_param*sizeof(CCTK_REAL));
    val_odd_Re =(CCTK_REAL *) malloc(*max_det_no_param*sizeof(CCTK_REAL));
    val_even_Im=(CCTK_REAL *) malloc(*max_det_no_param*sizeof(CCTK_REAL));
    val_odd_Im =(CCTK_REAL *) malloc(*max_det_no_param*sizeof(CCTK_REAL));

    for (il=*l_min; il<=*l_max; il=il+*l_step) {
      for (im=*m_min; im<=*m_max; im=im+*m_step) {
        if ((im<-il) || (im>+il)) continue;
        for (idet=0; idet<*max_det_no_param; idet++) {
          assert (idet>=0 && idet<*max_det_no_param);
          index = idet+(il-1)*(*max_det_no_param)+(-*m_min+im)*(*max_det_no_param)*(*l_max);
          assert (index>=0 && index<l_mode*(m_mode+1));

          val_even_Re[idet]=data_Qeven_Re[index];
          val_odd_Re[idet]=data_Qodd_Re[index];
          val_even_Im[idet]=data_Qeven_Im[index];
          val_odd_Im[idet]=data_Qodd_Im[index];
          if (verbose>5) {
            CCTK_VInfo(CCTK_THORNSTRING,
                  "idet: %d index: %d Re Qeven: %16.8g, Qodd: %16.8g",
                  (int)idet,(int)index,(double)val_even_Re[idet],(double)val_odd_Re[idet]);
            CCTK_VInfo(CCTK_THORNSTRING,
                  "idet: %d index: %d Im Qeven: %16.8g, Qodd: %16.8g",
                  (int)idet,(int)index,(double)val_even_Im[idet],(double)val_odd_Im[idet]);
          }
        }
        snprintf(string_var, sizeof string_var,
                 "Qeven_Re_l%d_m%d",(int)il,(int)im);
        WavExtrL_WriteScalar(cctkGH,vindex_Qeven_Re,
                            val_even_Re,maximum_detector_number,
                            string_var,detector_radius);
        snprintf(string_var, sizeof string_var,
                 "Qodd_Re_l%d_m%d",(int)il,(int)im);
        WavExtrL_WriteScalar(cctkGH,vindex_Qodd_Re,
                            val_odd_Re, maximum_detector_number,
                            string_var,detector_radius);
        if (verbose>3) CCTK_VInfo(CCTK_THORNSTRING,
                "   Detector 0: (l,m)=(%d,%d) Re Qeven: %16.8g Qodd: %16.8g",
                (int)il,(int)im,(double)val_even_Re[0],(double)val_odd_Re[0]);
        snprintf(string_var, sizeof string_var,
                 "Qeven_Im_l%d_m%d",(int)il,(int)im);
        WavExtrL_WriteScalar(cctkGH,vindex_Qeven_Im,
                            val_even_Im,maximum_detector_number,
                            string_var,detector_radius);
        snprintf(string_var, sizeof string_var,
                 "Qodd_Im_l%d_m%d",(int)il,(int)im);
        WavExtrL_WriteScalar(cctkGH,vindex_Qodd_Im,
                            val_odd_Im, maximum_detector_number,
                            string_var,detector_radius);
        if (verbose>3) CCTK_VInfo(CCTK_THORNSTRING,
                "   Detector 0: (l,m)=(%d,%d) Im Qeven: %16.8g Qodd: %16.8g",
                (int)il,(int)im,(double)val_even_Im[0],(double)val_odd_Im[0]);
      }
    }
    free(val_even_Re);
    free(val_odd_Re);
    free(val_even_Im);
    free(val_odd_Im);
  }

}



void WavExtrL_TimerInfo(CCTK_ARGUMENTS)
{
  cTimerData *info_interp, *info_schw, *info_MoncriefQ;
  CCTK_INT ierr;

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  if (*do_nothing==1)
    return;

  info_interp = CCTK_TimerCreateData();
  ierr = CCTK_Timer("Interpolation",info_interp);

  info_schw = CCTK_TimerCreateData();
  ierr = CCTK_Timer("Schwarzschild",info_schw);

  info_MoncriefQ = CCTK_TimerCreateData();
  ierr = CCTK_Timer("MoncriefQ",info_MoncriefQ);

  if (write_timer_info == 1)
  {
    CCTK_INFO("  Interpolation Timer");
    WavExtrL_WriteTimer(info_interp);

    CCTK_INFO("  Schwarzschild Timer (Integration 1)");
    WavExtrL_WriteTimer(info_schw);

    CCTK_INFO("  MoncriefQ Timer (Integration 2)");
    WavExtrL_WriteTimer(info_MoncriefQ);
  }

  /* cleanup */
  ierr = CCTK_TimerDestroyData(info_interp);
  ierr = CCTK_TimerReset("Interpolation");
  ierr = CCTK_TimerDestroyData(info_schw);
  ierr = CCTK_TimerReset("Schwarzschild");
  ierr = CCTK_TimerDestroyData(info_MoncriefQ);
  ierr = CCTK_TimerReset("MoncriefQ");
}

static void WavExtrL_WriteTimer(cTimerData *info)
{
  int i;

  for (i=0; i < info->n_vals; i++)
  {
    switch (info->vals[i].type)
    {
      case val_int:
        CCTK_VInfo(CCTK_THORNSTRING,"    %s %d %s",
        info->vals[i].heading, info->vals[i].val.i, info->vals[i].units);
        break;
      case val_long:
        CCTK_VInfo(CCTK_THORNSTRING,"    %s %ld %s",
        info->vals[i].heading, info->vals[i].val.l, info->vals[i].units);
        break;
      case val_double:
        CCTK_VInfo(CCTK_THORNSTRING,"    %s %.3f %s",
        info->vals[i].heading, (double)info->vals[i].val.d, info->vals[i].units);
        break;
      default:
        CCTK_WARN(1,"Unknown data type for timer info");
        break;
    }
  }
}

/* write scalars out to disk */
static int WavExtrL_WriteScalar(const cGH *GH,
                        CCTK_INT   vindex,
                        CCTK_REAL  *value,
                        CCTK_INT   max,
                        const char *alias,
                        const CCTK_REAL  *detector_radii)
{
  FILE *file;
  extractGH *myGH;  char *filename;
  const char *file_extension;
  char format_str_real[15], format_str_int[15], format_str_single[15];
  int i;

  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc (GH) != 0)
  {
    return (0);
  }

  /* set the output format string for the desired notation */
  snprintf (format_str_real, sizeof format_str_real,
            "%%%s\t%%%s\n", out_format, out_format);
  snprintf (format_str_int, sizeof format_str_int,
            "%%%s\t%%d\n", out_format);

  snprintf (format_str_single, sizeof format_str_single,
            "%%%s\t", out_format);

  /* get the GH extensions for Extract */
  myGH = (extractGH *) CCTK_GHExtension (GH, "WaveExtractL");

  /* set the output file extension according to the output style */
  file_extension = CCTK_Equals (out_style, "gnuplot") ? ".asc" : ".xg";

  filename = (char *) malloc (strlen (myGH->out_dir) + strlen (alias) +
                                strlen (file_extension) +1 );
  sprintf (filename, "%s%s%s", myGH->out_dir, alias, file_extension);

  if(max==-1)
  {
    /* create/reopen the file */
    file = WavExtrL_OpenFile(GH, vindex, filename, "tl", "Scalar value", 
                            -1, alias, NULL);
    if (file)
    {
      fprintf (file, format_str_real, GH->cctk_time, (double) *value);
      fclose (file);
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open scalar output file '%s'",
                  filename);
    }
  }
  else /* multiple detectors in one file */
  {
    if (CCTK_EQUALS(out_style, "xgraph"))
    {
      CCTK_ParameterSet(out_style,"WaveExtractL","gnuplot");
    }
    /* create/reopen the file */
    file = WavExtrL_OpenFile(GH, vindex, filename, "tl", "Scalar value",
                                max,alias,detector_radii);
    if (file)
    {
      fprintf (file, "\"Time = %g\n",(double)GH->cctk_time);
      /* this is slow, but we don't care */
      for (i=0; i<max; i++) {
        fprintf (file, format_str_real, (double) detector_radii[i],
                                        (double) value[i]);
      }
      fprintf (file,"\n\n");
      fclose (file);
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open scalar output file '%s'",
                  filename);
    }
  }

  /* clean up */
  free (filename);
  return (0);
}


static FILE *WavExtrL_OpenFile (const cGH *GH,
                               int vindex,
                               const char *filename,
                               const char *slicename,
                               const char *description,
                               CCTK_INT   number_sets,
                               const char *aliasname,
                               const CCTK_REAL  *detector_radii)
{
  int first_time_through;
  FILE *file;
  extractGH *myGH;
  const ioGH *ioUtilGH;
  char *fullname;
  char comment_char, buffer[10000];
  ioAdvertisedFileDesc advertised_file;
  DECLARE_CCTK_PARAMETERS

  /* get the GH extension handles for IOUtil and Extract */
  myGH = (extractGH *) CCTK_GHExtension (GH, "WaveExtractL");
  ioUtilGH = (const ioGH *) CCTK_GHExtension (GH, "IO");

  /* get the variable's full name */
  fullname = CCTK_FullName (vindex);

  /* create the output file the first time through
     If restart from recovery, existing files are opened in append mode. */
  first_time_through = GetNamedData (myGH->filenameList, filename)==NULL;
  file = fopen (filename,
                ioUtilGH->recovered || ! first_time_through ?  "a" : "w");
  if (! file)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "WavExtrL_OpenFile: Cannot open output file '%s'", filename);
  }
  else if (first_time_through)
  {
    if (CCTK_Equals (out_style, "gnuplot"))
    {
      comment_char = '#';
      advertised_file.mimetype = "application/gnuplot";
    }
    else
    {
      comment_char = '"';    /* this is for xgraph */
      advertised_file.mimetype = "application/x-graph";
    }

    /* just store a non-NULL pointer in database */
    StoreNamedData (&myGH->filenameList, filename, (void *) 1);

    /* advertise the file for downloading */
    advertised_file.slice = slicename;
    advertised_file.thorn = CCTK_THORNSTRING;
    advertised_file.varname = fullname;
    advertised_file.description = description;

    IOUtil_AdvertiseFile (GH, filename, &advertised_file);

    /* don't print the header again after recovery */
    if (! ioUtilGH->recovered)
    {
      /* write the file info and the header */
      if (CCTK_Equals (out_fileinfo, "parameter filename") ||
          CCTK_Equals (out_fileinfo, "all"))
      {
        buffer[0] = 0;
        CCTK_ParameterFilename (sizeof (buffer), buffer);
        fprintf (file, "%cParameter file %s\n", comment_char, buffer);
      }
      if (CCTK_Equals (out_fileinfo, "creation date") ||
          CCTK_Equals (out_fileinfo, "all"))
      {
        buffer[0] = 0;
        Util_CurrentDate (sizeof (buffer), buffer);
        fprintf (file, "%cCreated %s ", comment_char, buffer);
        Util_CurrentTime (sizeof (buffer), buffer);
        fprintf (file, "%s\n", buffer);
      }
      if (number_sets==-1)
      {
        if (CCTK_Equals (out_fileinfo, "axis labels") ||
            CCTK_Equals (out_fileinfo, "all"))
        {
          fprintf (file, "%cx-label time\n", comment_char);
          fprintf (file, "%cy-label %s\n", comment_char, advertised_file.varname);
        }
      }
      else
      {
        if (CCTK_Equals (out_fileinfo, "axis labels") ||
            CCTK_Equals (out_fileinfo, "all"))
        {
          fprintf (file, "%cx-label x\n", comment_char);
          fprintf (file, "%cy-label %s\n", comment_char,advertised_file.varname);
          fprintf(file, "\n\n");
        }
      }
    }
  }

  free (fullname);

  return (file);
}


 /*@@
   @routine   WavExtrL_SetupGH
   @date      Sat Feb 6 1999
   @author    Gabrielle Allen
   @desc
              Allocates and sets up Extract's GH extension structure
   @enddesc

   @calls     CCTK_RegisterIOMethod
              CCTK_RegisterIOMethodOutputGH
              CCTK_RegisterIOMethodTimeToOutput
              CCTK_RegisterIOMethodTriggerOutput

   @var       config
   @vdesc     the CCTK configuration as provided by the flesh
   @vtype     tFleshConfig *
   @vio       unused
   @endvar
   @var       conv_level
   @vdesc     the convergence level
   @vtype     int
   @vio       unused
   @endvar
   @var       GH
   @vdesc     Pointer to CCTK grid hierarchy
   @vtype     cGH *
   @vio       in
   @endvar

   @returntype  void *
   @returndesc
                pointer to the allocated GH extension structure
   @endreturndesc
@@*/
static void *WavExtrL_SetupGH (tFleshConfig *config, int conv_level, cGH *GH)
{
  int i;
  extractGH *myGH;
  const char *my_out_dir;
  DECLARE_CCTK_PARAMETERS

  CCTK_INFO("setup GH");

  /* suppress compiler warnings about unused parameters */
  (void) (config + 0);
  (void) (conv_level + 0);
  (void) (GH + 0);

  /* allocate the GH extension and its components */
  myGH = (extractGH *) malloc (sizeof (extractGH));
  if (myGH)
  {
    i = CCTK_RegisterIOMethod ("WaveExtractL");
    CCTK_RegisterIOMethodTimeToOutput (i, WavExtrL_TimeForOutput);

    i = CCTK_NumVars ();

    myGH->filenameList = NULL;
    myGH->this_time = 0.0;

    /* get the name of Extract's output directory */
    my_out_dir = out_dir;
    if (*my_out_dir == 0)
    {
      my_out_dir = io_out_dir;
    }

    /* skip the directory pathname if output goes into current directory */
    if (strcmp (my_out_dir, "."))
    {
      i = strlen (my_out_dir);
      myGH->out_dir = (char *) malloc (i + 2);
      strcpy (myGH->out_dir, my_out_dir);
      myGH->out_dir[i] = '/';
      myGH->out_dir[i+1] = 0;
    }
    else
    {
      myGH->out_dir = strdup ("");
    }

    /* create the output dir */
    if (*myGH->out_dir && CCTK_MyProc (GH) == 0)
    {
      i = IOUtil_CreateDirectory (GH, myGH->out_dir, 0, 0);
      if (i < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Couldn't create Extract output directory "
                    "'%s'", myGH->out_dir);
      }
      else if (i >= 0 && verbose)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Output directory is '%s'",
                    myGH->out_dir);
      }
    }

    CCTK_VInfo(CCTK_THORNSTRING,"Setting out_every to %d",(int)out_every);
    myGH->out_every=out_every;
  }
  else {
    CCTK_WARN(1,"unable to get GH");
  }

  return (myGH);
}



/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine    WavExtrL_CheckSteerableParameters
   @date       Fri Dec 20 2002
   @author     Gabrielle Allen
   @desc
               Re-evaluates 'Extract::out_every' and/or 'IO::out_every'
               resp. to set myGH->out_every to the frequency of Extract
               output. Copied from Thomas Radke routines in other IO thorns.
   @enddesc
   @calls      CCTK_ParameterQueryTimesSet

   @var        myGH
   @vdesc      Pointer to IOBasic's GH extension
   @vtype      iobasicGH *
   @vio        in
   @endvar
@@*/
void WavExtrL_CheckSteerableParameters (extractGH *myGH)
{
  int out_old;
  int i;
  char parname[256];
  char out_everyname[256];
  DECLARE_CCTK_PARAMETERS

  if (verbose>2) CCTK_INFO("Checking steerable out_every parameter");

  /* how often to output */

  out_old = myGH->out_every;
  myGH->out_every = out_every >= 0 ? out_every : io_out_every;

  /* report if frequency changed */
  if (myGH->out_every != out_old )
  {
    if (myGH->out_every > 0)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Periodic output changed"
                  "from %d to every %d iterations", out_old, myGH->out_every);
      /* change all individual detectors */
      snprintf(out_everyname,sizeof out_everyname,"%d",(int)out_every);
      for (i=0;i<=maximum_detector_number;i++) {
        /* can't set array parameters in one go, build up name */
        snprintf(parname,sizeof parname,"out_every_det[%d]",i);
        if (verbose>3) printf("setting %s\n",parname);
        CCTK_ParameterSet(parname, CCTK_THORNSTRING, out_everyname);
      }
    }
    else
    {
      CCTK_INFO ("Periodic output turned off");
    }
  }

  /* return if there's nothing to do */
  if (myGH->out_every <= 0)
  {
    return;
  }

}

