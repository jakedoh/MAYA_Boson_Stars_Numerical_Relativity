 /*@@
   @header    ioHDF5GH.h
   @date      Jun 20 2000
   @author    Thomas Radke
   @desc
              The extensions to the GH structure from IOHDF5.
   @enddesc
   @version   $Header$
 @@*/

#ifndef _IOHDF5_IOHDF5GH_H_
#define _IOHDF5_IOHDF5GH_H_ 1

#include "StoreNamedData.h"
#include "CactusPUGHIO/IOHDF5Util/src/ioHDF5UtilGH.h"

/* IOHDF5 GH extension structure */
typedef struct
{
  /* default number of times to output */
  int out_every_default;

  /* number of times to output for each variable */
  CCTK_INT *out_every;

  /* the last iteration output for each variable */
  int *out_last;

  /* list of variables to output */
  char *out_vars;

  /* I/O request description list (for all variables) */
  ioRequest **requests;

  /* directory in which to output */
  char *out_dir;

  /* filename database for opened files */
  pNamedData *open_output_files;

  /* timer array for checkpointing/recovery */
  int timers[IOHDF5_NUM_TIMERS];

  /* flag to indicate request for timer output */
  int print_timing_info;

  /* ring buffer for list of successfully created cp files */
  int    checkpoint_keep;
  int    cp_filename_index;
  char **cp_filename_list;

  /* iteration number of the last checkpoint */
  int last_checkpoint_iteration;

  /* stop on I/O parameter parsing errors ? */
  int stop_on_parse_errors;

} ioHDF5GH;

#ifdef __cplusplus
extern "C"
{
#endif

/* prototypes of functions to be registered as IOHDF5's IO method */
int IOHDF5_OutputGH (const cGH *GH);
int IOHDF5_TriggerOutput (const cGH *GH, int);
int IOHDF5_TimeFor (const cGH *GH, int);
int IOHDF5_OutputVarAs (const cGH *GH, const char *var, const char *alias);
int IOHDF5_Recover (cGH *GH, const char *basefilename, int called_from);

/* other function prototypes */
int IOHDF5_Write (const cGH *GH, int vindex, const char *alias);
int IOHDF5_WriteIsosurface (const cGH *GH,
                            const char *filename,
                            const char *GVname,
                            CCTK_INT iteration,
                            CCTK_INT timelevel,
                            CCTK_REAL isoval,
                            CCTK_REAL minval,
                            CCTK_REAL maxval,
                            int nTriangles,
                            const CCTK_INT4 *triangles,
                            int nVertices,
                            const CCTK_REAL4 *vertices);
void IOHDF5_CheckSteerableParameters (const cGH *GH, ioHDF5GH *myGH);
#ifdef __cplusplus
} // extern "C"
#endif

#endif  /* _IOHDF5_IOHDF5GH_H_ */
