
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
   @file      ResetCurrDet.c
   @date      16 Jul 2003
   @author    Frank Herrmann
   @desc 
              Reset the current_detector value, this is needed for the while loop in
              the scheduler.
   @enddesc 
   @version $Id: ResetCurrDet.c 37 2008-02-14 03:30:05Z schnetter $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void WavExtrL_ResetCurrDet(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtrL_Reset_CurrDet
   @date       16 Jul 2003
   @author     Frank Herrmann
   @desc 
               Reset the current_detector value, this is needed for the while loop in
               the scheduler.
   @enddesc 
   @calls     
@@*/
void WavExtrL_ResetCurrDet(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>0)
    CCTK_VInfo(CCTK_THORNSTRING,"Calling WaveExtractL at time: %f",(double)cctkGH->cctk_time);

  if (verbose>2)
    CCTK_INFO("Reset the current_detector value");

  if (verbose >3)
    CCTK_VInfo(CCTK_THORNSTRING,"iteration %d detector %d out_every_det[] %d",cctk_iteration,(int)*current_detector,(int)my_out_every_det[*current_detector]);

  if (cctk_iteration != 0) {
    if (cctk_iteration % my_out_every_det[*current_detector] != 0)
    {
      *do_nothing=1;
    }
    else if (Cauchy_time_ID!=0) {
      if((cctk_iteration-1) % my_out_every_det[*current_detector] == 0) {
        *do_nothing=0;
      }
      else if ((cctk_iteration+1) % my_out_every_det[*current_detector] == 0) {
        *do_nothing=0;
      }
    }

    if (cctk_iteration % my_out_every_det[*current_detector] == 0) {
      *do_nothing=0;
    }
  }
  else
    *do_nothing = 0;

  if (*do_nothing == 1) {
    return;
  }

  *current_detector=maximum_detector_number;

}
