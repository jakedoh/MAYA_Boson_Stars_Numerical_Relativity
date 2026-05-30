 /*@@
   @header    extractGH.h
   @date      Friday 18th September 1999
   @author    Gabrielle Allen
   @desc
              The extensions to the GH structure from WaveExtract.
   @enddesc
   @version $Header$
 @@*/

#include "StoreNamedData.h"

typedef struct EXTRACT_GH
{
  CCTK_REAL this_time;
  CCTK_REAL last_time; 

  /* how often to output */
  int  out_every;

  /* directory in which to place scalar output */
  char *out_dir;

  /* database for names of output files that were already created */
  pNamedData *filenameList;

} extractGH;


/* prototypes of functions to be registered */
int WavExtr_TimeForOutput (const cGH *GH, int vindex);
void WavExtr_CheckSteerableParameters (extractGH *myGH);
