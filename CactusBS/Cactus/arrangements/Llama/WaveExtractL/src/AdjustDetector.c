
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
   @file      AdjustDetector.c
   @date      Mon Nov 25 13:01:43 2002
   @author    Frank Herrmann
   @desc 
              Adjust next detector, increase the variable current_detector
   @enddesc 
   @version  $Id: AdjustDetector.c 6 2004-06-30 18:08:24Z herrmann $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void WavExtrL_AdjustDetector(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtrL_AdjustDetector
   @date       Mon Nov 25 13:02:44 2002
   @author     Frank Herrmann
   @desc 
               Adjust next detector. This function decreases the variable 
               current_detector by one.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
void WavExtrL_AdjustDetector(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>4)
    CCTK_INFO("Adjust next detector");

  if (*do_nothing == 1)
    *do_nothing = 0;

  if (*current_detector != 0) 
    *current_detector=*current_detector-1;

  if (*current_detector < 0)
    *current_detector=0;

  if (verbose >2)
    CCTK_VInfo(CCTK_THORNSTRING,"Done with detector %d, next will go to Detector No. %d",(int)*current_detector +1,(int)*current_detector);
}
