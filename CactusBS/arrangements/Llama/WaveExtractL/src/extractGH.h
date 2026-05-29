
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
   @header    extractGH.h
   @date      Friday 18th September 1999
   @author    Gabrielle Allen
   @desc
              The extensions to the GH structure from WaveExtract.
   @enddesc
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
int WavExtrL_TimeForOutput (const cGH *GH, int vindex);
void WavExtrL_CheckSteerableParameters (extractGH *myGH);
