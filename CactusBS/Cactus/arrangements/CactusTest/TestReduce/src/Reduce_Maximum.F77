#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine TestReduceF_Maximum(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS

      integer handle
      integer index
      integer ierr
      CCTK_REAL reduction_value
      CCTK_INT ireduction_value
      
      call CCTK_INFO("TESTING MAXIMUM REDUCTIONS FROM FORTRAN")

c     3D GF

      call CCTK_VarIndex(index,"TestReduce::phi_gf3")
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
	 call CCTK_Reduce (ierr,cctkGH, 
     &        0,handle, 1, 
     &        CCTK_VARIABLE_REAL,
     &        CCTK_PointerTo (reduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_gf3 is ",reduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if
      
      call CCTK_VarIndex(index,"TestReduce::phi_igf3")
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
	 call CCTK_Reduce (ierr,cctkGH, 
     &        0,handle, 1, 
     &        CCTK_VARIABLE_INT,
     &        CCTK_PointerTo (ireduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_igf3 is ",ireduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if
      
c     2D GF 

      call CCTK_VarIndex(index,"TestReduce::phi_gf2")      
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
         call CCTK_Reduce (ierr,cctkGH,0, 
     &        handle,1, 
     &        CCTK_VARIABLE_REAL,
     &        CCTK_PointerTo (reduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_gf2 is",reduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if

      call CCTK_VarIndex(index,"TestReduce::phi_igf2")      
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
         call CCTK_Reduce (ierr,cctkGH,0, 
     &        handle,1, 
     &        CCTK_VARIABLE_INT,
     &        CCTK_PointerTo (ireduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_igf2 is",ireduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if


c     1D GF 

      call CCTK_VarIndex(index,"TestReduce::phi_gf1")
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
         call CCTK_Reduce (ierr,cctkGH,0, 
     &        handle,1, 
     &        CCTK_VARIABLE_REAL,
     &        CCTK_PointerTo (reduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_gf1 is",reduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if

      call CCTK_VarIndex(index,"TestReduce::phi_igf1")
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
         call CCTK_Reduce (ierr,cctkGH,0, 
     &        handle,1, 
     &        CCTK_VARIABLE_INT,
     &        CCTK_PointerTo (ireduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH) .eq. 0) then
            print *,"Maximum value of phi_igf1 is",ireduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if

c     Grid Scalar

      call CCTK_VarIndex(index,"TestReduce::phi_scalar")  
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
         call CCTK_Reduce (ierr,cctkGH, 
     &        0,handle,1, 
     &        CCTK_VARIABLE_REAL,
     &        CCTK_PointerTo (reduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_scalar is ",reduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if

      call CCTK_VarIndex(index,"TestReduce::phi_iscalar")  
      call CCTK_ReductionHandle(handle,"maximum")

      if (handle .gt. 0) then
         call CCTK_Reduce (ierr,cctkGH, 
     &        0,handle,1, 
     &        CCTK_VARIABLE_INT,
     &        CCTK_PointerTo (ireduction_value),
     &        1,index)
         if (ierr .eq. 0 .and. CCTK_MyProc(cctkGH).eq.0) then
            print *,"Maximum value of phi_scalar is ",ireduction_value
         end if
      else
         call CCTK_WARN (1,"TestReduce: Invalid reduction operator")
      end if

      return
      end

