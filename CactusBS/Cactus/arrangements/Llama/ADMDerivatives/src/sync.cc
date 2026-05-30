
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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


extern "C" void ADMDerivatives_radial_SelectBC(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   CCTK_INT faces;
  
   int ierr = 0;

   faces = CCTK_ALL_FACES;

   if (store_radial_derivatives)
   {
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, 1, -1,
                                           "ADMDerivatives::dr_lapse", "None");
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, 1, -1,
                                           "ADMDerivatives::dr_shift", "None");
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, 1, -1,
                                           "ADMDerivatives::dr_metric", "None");
   }
   
   if (ierr)
      CCTK_WARN(1, "Error selecting BCs!");
}



extern "C" void ADMDerivatives_cartesian_SelectBC(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   CCTK_INT faces;
  
   int ierr = 0;

   faces = CCTK_ALL_FACES;

   if (store_cartesian_derivatives)
   {
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, 1, -1,
                                           "ADMDerivatives::dx_vars", "None");
   }
   
   if (ierr)
      CCTK_WARN(1, "Error selecting BCs!");
}



extern "C" void ADMDerivatives_time_SelectBC(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   CCTK_INT faces;
  
   int ierr = 0;

   faces = CCTK_ALL_FACES;

   if (store_time_derivatives)
   {
      ierr += Boundary_SelectGroupForBC(cctkGH, faces, 1, -1,
                                           "ADMDerivatives::dt_metric", "None");
   }
   
   if (ierr)
      CCTK_WARN(1, "Error selecting BCs!");
}

