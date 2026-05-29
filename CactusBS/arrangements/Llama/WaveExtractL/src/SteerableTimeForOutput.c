
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
  @file      Output.c
  @date      April 11th 2002
  @author    Frank Herrmann
  @desc
             Check if it is the right time for WaveExtract.
             
  @enddesc
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"

#include "extractGH.h"


static int WaveExtract_ncall=0; /* number of calls to WaveExtract */

void WavExtrL_CheckSteerableParameters(extractGH *myGH);


/* prototype */
int WavExtrL_TimeForOutput(const cGH *GH, int vindex);

/*@@
  @routine    WavExtrL_TimeForOutput
  @date       April 11th 2002
  @author     Frank Herrmann
  @desc
              Check if it is time for output
  @enddesc
  @calls      CheckSteerableParameters

  @var        GH
  @vdesc      Pointer to CCTK GH
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        vindex
  @vdesc      index of variable to check for output
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              true/false (1 or 0) if analysis should be called
  @endreturndesc
@@*/
int WavExtrL_TimeForOutput(const cGH *GH, int vindex)
{
  extractGH *myGH;
  int retval=1;            /* by default we do analysis, next we check when not */
  int WaveExtract_after;  /* after which iteration we start analysis */

  DECLARE_CCTK_PARAMETERS


  /* FIXME : There is a lot of crap in this file and in IO.c as well */


  if (verbose > 6)
    CCTK_INFO("Triggering");

  /* FIXME : STEERABLE PARAMETERS */
  /* get myGH extract GH */
  myGH=(extractGH *) CCTK_GHExtension (GH, "WaveExtractL");
  WavExtrL_CheckSteerableParameters (myGH);

  /* check we should do output at all */
  if(out_every == 0) { 
    if (verbose > 5)
      CCTK_INFO("You set out_every=0");
    retval=0;
  }
  else
  {
    /* check we have the right variable */
    if (vindex==CCTK_VarIndex("WaveExtractL::triggervar"))
    {
      if (verbose > 1 )
        CCTK_INFO("Starting WaveExtract");

      /* sanity checks */
      if(start_iteration <0 && start_time <0)
      {
        CCTK_WARN(1,
          "you need to specify start_time or start_iteration");
        retval=0;
        return retval;
      }

      if(start_iteration >= 0 && start_time >= 0)
      {
        CCTK_WARN(1,
          "specify only one of start_time or start_iteration");
        retval=0;
        return retval;
      }

      /* have we reached the start iteration? */
      if (start_iteration >= 0) 
      {
        WaveExtract_after=start_iteration;
        if (GH->cctk_iteration <WaveExtract_after) 
          retval=0;
      }

      /* have we reached the start time? */
      if (start_time >= 0)
      {
        if (GH->cctk_time <start_time) 
          retval=0;

        if (WaveExtract_ncall == 0) 
          WaveExtract_after=GH->cctk_iteration;
      }
 
      /* are we at the right iteration? */ 
      if ( (GH->cctk_iteration-WaveExtract_after) %out_every) 
        retval=0;

      if (retval==1) 
      {
        WaveExtract_ncall=0;  /* store analysis in extract_ncall */
        /*myGH->last_time = myGH->this_time;
        myGH->this_time = GH->cctk_time; */
      }
    }
    /* not the right variable */
    else
    {
      if (verbose > 9 && retval !=0 )
        CCTK_INFO("triggered with the wrong variable");
      retval = 0;
    }
  }

  return retval;
}

