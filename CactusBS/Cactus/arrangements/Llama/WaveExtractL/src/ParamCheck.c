
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
   @file      ParamCheck.c
   @date      April 26 2002
   @author    Gabrielle Allen
   @desc 
              Check parameters for WaveExtract
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>

void WavExtrL_ParamCheck(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtrL_ParamCheck
   @date       April 26 2002
   @author     Gabrielle Allen
   @desc 
               Check parameters for WaveExtract
   @enddesc 
   @calls     
@@*/
void WavExtrL_ParamCheck(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>3)
    CCTK_INFO("Checking parameters");

  if(! (CCTK_EQUALS(metric_type, "physical") || 
        CCTK_EQUALS(metric_type, "static conformal")))
  {
    CCTK_PARAMWARN("WaveExtract only works currently with metric_type \"static conformal\" or \"physical\"");
  }
  
  int store_radial_derivs = *(CCTK_INT const*) CCTK_ParameterGet("store_radial_derivatives", "ADMDerivatives", NULL);
  if (!store_radial_derivs)
     CCTK_PARAMWARN("Radial derivatives need to be calculated in ADMDerivatives!");

  if (origin_x != 0.0 || origin_y != 0.0 || origin_z !=0)
     CCTK_PARAMWARN("Radial derivatives as calculated in ADMDerivatives assume origin(x,y,z) = 0.0!");

  if (!check_rmax && (Cauchy_radius_start_factor > 0. || Cauchy_radius_end_factor > 0))
     CCTK_PARAMWARN("When using Cauchy_radius_start_factor or Cauchy_radius_end_factor, check_rmax must also be set");
}
