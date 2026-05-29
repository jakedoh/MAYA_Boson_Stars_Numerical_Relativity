 /*@@
   @header    ioJpegGH.h
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              The extensions to the GH structure from IOJpeg.
   @enddesc
   @version   $Header$
 @@*/

#ifndef IOJPEG_IOJPEGGH_H_
#define IOJPEG_IOJPEGGH_H_ 1

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct IOJpegGH
{
  /* default number of times to output */
  int out_every_default;

  /* number of times to output every variable */
  CCTK_INT *out_every;

  /* the last iteration output for every variable */
  int *out_last;

  /* list of variables to output */
  char *out_vars;

  /* directories in which to output */ 
  char *out_dir;

  /* for 2D planes, we define the index where to locate the plane
     sp2xyz[maxdim][perpendicular direction] */
  int **sp2xyz;

  /* stop on I/O parameter parsing errors ? */
  int stop_on_parse_errors;

} ioJpegGH;


/* prototypes of functions to be registered as I/O method */
int IOJpeg_OutputGH (const cGH *GH);
int IOJpeg_OutputVarAs (const cGH *GH, const char *fullname, const char *alias);
int IOJpeg_TimeFor (const cGH *GH, int vindex);
int IOJpeg_TriggerOutput (const cGH *GH, int vindex);

/* other function prototypes */
int IOJpeg_Write (const cGH *GH, CCTK_INT vindex, const char *alias);
void IOJpeg_CheckSteerableParameters (ioJpegGH *myGH);

/* routines called from JPEG.c */
int WriteJPEGToFileRGB (int nx, int ny, void *data, int Quality, FILE* outfile);
int WriteJPEGToMemoryRGB (int nx, int ny, void *data, int Quality, char *buffer,
                          int buffersize);
void AutoColorDataSlice (int nx, int ny, const CCTK_REAL *datain,
                         unsigned char *dataout, CCTK_REAL min, CCTK_REAL max,
                         CCTK_REAL bias, int rdfac);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* IOJPEG_IOJPEGGH_H_ */
