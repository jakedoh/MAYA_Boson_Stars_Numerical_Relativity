#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_IOMethods.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
inline void CREATE_OUTDIR(CCTK_ARGUMENTS,const char *input_dir,char *&actual_dir)
  {
    DECLARE_CCTK_PARAMETERS;
    const char *_dir = input_dir;


    /* check whether "dir" was set; if not default to "IO::out_dir" */
    if (*_dir == 0)
      {
        _dir = out_dir;
      }

    /* omit the directory name if it's the current working dir */
    if (strcmp (_dir, ".") == 0)
      {
        actual_dir = strdup ("");
      }
    else
      {
        actual_dir = (char *)malloc (strlen (_dir) + 2);
        sprintf (actual_dir, "%s/", _dir);
      }

    /* create the directory */
    int izzle = IOUtil_CreateDirectory (cctkGH, actual_dir, 0, 0);
    if (izzle < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "VolumeIntegrals_GRMHD:file_output_routines.C Problem creating directory '%s' "
                    "for output", actual_dir);
      }
    else if (izzle >= 0 && verbose==1)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Outputting to directory '%s'",
                    actual_dir);
      }
  }

inline void OUTDIR_NAME(CCTK_ARGUMENTS,const char *input_dir,char *&actual_dir)
{
    DECLARE_CCTK_PARAMETERS;
    const char *_dir = input_dir;


    /* check whether "dir" was set; if not default to "IO::out_dir" */
    if (*_dir == 0)
      {
        _dir = out_dir;
      }

    /* omit the directory name if it's the current working dir */
    if (strcmp (_dir, ".") == 0)
      {
        actual_dir = strdup ("");
      }
    else
      {
        actual_dir = (char *)malloc (strlen (_dir) + 2);
        sprintf (actual_dir, "%s/", _dir);
      }
}

void particle_tracerET_file_output_routine_Startup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  char *actual_dir;
  CREATE_OUTDIR(CCTK_PASS_CTOC,outdir,actual_dir);
}

void particle_tracerET_file_output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if(update_RK4_freq==0) return;

  /* Only output on the zeroth MPI process, after a full RK4 step or at cctk_iteration==0.
   *    Note that each RK4 step first goes from t=0 to t=dt/2, then completes at t=dt.
   *    Thus, a full RK4 step is completed every 2*update_RK4_freq. */
  if( CCTK_MyProc(cctkGH)==0 && (cctk_iteration==start_tracing_particles_iteration || cctk_iteration%(2*update_RK4_freq)==0) ) {
    char *actual_dir;
    OUTDIR_NAME(CCTK_PASS_CTOC,outdir,actual_dir);
  
    /* Next loop over all integrals, output a file for each */

    /* allocate filename string buffer and build the filename */
    char *filename = (char *)malloc (strlen (actual_dir) +
				     strlen ("particles.asc") + 10); // <- Tiny amount of memory.
    sprintf (filename, "%sparticles.asc", actual_dir);
    FILE *file = fopen (filename,"a+");
    if (! file) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "particle_tracerET: Cannot open output file '%s'", filename);
      exit(1);
    } else {

      fseek(file, 0, SEEK_END);
      long size = ftell(file);
      /* First print header if file size is zero. */
      if(size==0) {
	char *header_buffer = (char *) malloc(sizeof(char)*500); // <- Tiny amount of memory.
	sprintf(header_buffer,"# Col. 1: Time (cctk_time). Subsequent columns: x,y,z coordinate of particle i\n");
	fprintf(file,"%s",header_buffer);
        free(header_buffer);
      }
      /* Next print one line of data, corresponding to a single point in time. */
      char *buffer = (char *) malloc(sizeof(char)*5000000); // <- Allocate 5MB; sizeof(char) is 1.
      sprintf(buffer,"%e ",cctk_time);
      for(int which_particle=0;which_particle<num_particles;which_particle++) {
	if(which_particle!=num_particles-1) {
          sprintf(buffer, "%s%e %e %e ",buffer,particle_position_x[which_particle],particle_position_y[which_particle],particle_position_z[which_particle]);
        } else {
          sprintf(buffer, "%s%e %e %e", buffer,particle_position_x[which_particle],particle_position_y[which_particle],particle_position_z[which_particle]);
        }
      }
      fprintf(file,"%s\n",buffer);
      fclose(file);
      free(buffer);
    }
    free(filename);
  }
}
