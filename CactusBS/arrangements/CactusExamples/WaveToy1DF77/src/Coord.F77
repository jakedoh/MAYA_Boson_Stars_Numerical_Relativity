 /*@@
   @file      Coord.F77
   @date      June 13 2000
   @author    Gabrielle Allen
   @desc 
              Coordinates for the 1d wave equation solver
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


 /*@@
   @routine    WaveToy1DF77_RegisterCoords
   @date       June 13 2000
   @author     Gabrielle Allen 
   @desc 
     Routine registers the coordinates for WaveToy1Df77
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

      SUBROUTINE WaveToy1DF77_RegisterCoords(CCTK_ARGUMENTS)

      implicit none
  
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_INT, parameter :: izero = 0
      integer, parameter :: ik = kind(izero)
      
      integer retval,ierr
      CCTK_INT coord_system_handle

      retval = 0

      ierr= Coord_SystemRegister(cctkGH, 1_ik, "cart1d");

C      call CCTK_CoordRegisterSystem(ierr,1,"cart1d")

c      call CCTK_INFO("Getting System Handle")

      coord_system_handle = Coord_SystemHandle(cctkGH, "cart1d")

c      call CCTK_INFO("Setting Coord Handle")

      ierr =  Coord_CoordRegister(cctkGH, coord_system_handle, 1_ik, "x")

c      call CCTK_INFO("Setting Coord Type")

      call Util_TableSetString(ierr, coord_system_handle, "uniform", "TYPE")

C      call CCTK_CoordRegisterData(ierr,1,"wavetoy1df77::x1d","x","cart1d")

      if (ierr .lt. 0) then
         call CCTK_WARN(1,"Problem with registering coordinate x")
         retval = -1
      end if

c      call CCTK_INFO("Setting Default System")

      ierr = Coord_SetDefaultSystem(cctkGH, "cart1d");

      return
      end


 /*@@
   @routine    WaveToy1DF77_Coord
   @date       June 13 2000
   @author     Gabrielle Allen
   @desc 
               Coordinates for the 1d wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToy1DF77_Coord(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS

      INTEGER   i,ierr
      CCTK_REAL xmin,xmax,x_origin,lowerx,upperx
      CCTK_INT coord_handle
      INTEGER varindex

      xmin=0
      xmax=1

c     Grid spacing on coarsest grid 
      cctk_delta_space(1) = (xmax - xmin) / max(cctk_gsh(1)- 1, 1)
      x_origin = xmin

      do i=1,cctk_lsh(1)
         x(i) = cctk_delta_space(1)*(i-1+cctk_lbnd(1)) + x_origin
      end do

      lowerx = x_origin
      upperx = x_origin+cctk_delta_space(1)*(cctk_gsh(1)-1)

c      call CCTK_INFO("Getting Coord Handle")

      coord_handle = Coord_CoordHandle(cctkGH, "x", "cart1d")
      if (coord_handle .lt. 0) then
         call CCTK_WARN(0, "Error retrieving coordinate handle for x of cart1d")
      endif

c      call CCTK_INFO("Getting Variable Index")

      call CCTK_VarIndex(varindex, "wavetoy1df77::x")

c      call CCTK_INFO("Setting Physical Minimum")

      call Util_TableSetReal   (ierr,coord_handle, lowerx, "PHYSICALMIN")

c      call CCTK_INFO("Setting Computational Minimum")

      call Util_TableSetReal   (ierr,coord_handle, lowerx, "COMPMIN")

c      call CCTK_INFO("Setting Physical Maximum")

      call Util_TableSetReal   (ierr,coord_handle, upperx, "PHYSICALMAX")

c      call CCTK_INFO("Setting Computational Maximum")

      call Util_TableSetReal   (ierr,coord_handle, upperx, "COMPMAX")

c      call CCTK_INFO("Setting Coordinate Type")

      call Util_TableSetString (ierr,coord_handle, "uniform", "TYPE")

c      call CCTK_INFO("Setting Time dependency")

      call Util_TableSetString (ierr,coord_handle, "no", "TIMEDEPENDENT")

c      call CCTK_INFO("Setting datatype")

      call Util_TableSetString (ierr,coord_handle, "CCTK_REAL", "DATATYPE")

c      call CCTK_INFO("Setting index")

      call Util_TableSetInt    (ierr,coord_handle, varindex, "GAINDEX")

c      call CCTK_INFO("Setting delta")

      call Util_TableSetReal   (ierr,coord_handle, cctk_delta_space(1), "DELTA")

c      call CCTK_INFO("Finished Setting up Coordinates")

      return
      end 


