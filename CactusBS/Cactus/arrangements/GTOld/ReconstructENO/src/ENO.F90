 /*@@
   @file      ENO.F90
   @date      Sat Apr  6 17:37:56 2002
   @author    Ian Hawke
   @desc 
   Routines to set up the coefficient array and to perform one dimensional 
   ENO reconstruction of arbitrary order.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    ENO_Setup
   @date       Sat Apr  6 17:42:13 2002
   @author     Ian Hawke
   @desc 
   Sets up the coefficient array for ENO reconstruction. 
   Uses the notation of Shu, equation (2.21), in 
   ''High Order ENO and WENO Schemes for CFD''
   (see ThornGuide for full reference).
   One exception: (Shu) r -> (Here) i: avoiding name clash.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine ENO_Setup(CCTK_ARGUMENTS)

  USE ENO_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k, l, m, q, allocstat
  CCTK_REAL :: denominator, numerator, numerator_product

  if(.not.coeffs_allocated) then
     allocate(eno_coeffs(-1:eno_order - 1, 0:eno_order - 1), STAT=allocstat)

     if (allocstat .ne. 0) call CCTK_WARN(0, "Failed to allocate ENO coefficient arrays!")
     coeffs_allocated = .true.
  endif

  do i = -1, eno_order - 1
    do j = 0, eno_order - 1

      eno_coeffs(i, j) = 0.d0

      do m = j+1, eno_order

        denominator = 1.d0
        do l = 0, m-1
          denominator = denominator * dble(m - l)
        end do
        do l = m+1, eno_order
          denominator = denominator * dble(m - l)
        end do

        numerator = 0.d0
        do l = 0, m-1
          numerator_product = 1.d0
          do q = 0, l-1
            numerator_product = numerator_product * dble(i - q + 1)
          end do
          do q = l+1, m-1
            numerator_product = numerator_product * dble(i - q + 1)
          end do
          do q = m+1, eno_order
            numerator_product = numerator_product * dble(i - q + 1)
          end do
          numerator = numerator + numerator_product
        end do

        do l = m+1, eno_order
          numerator_product = 1.d0
          do q = 0, m-1
            numerator_product = numerator_product * dble(i - q + 1)
          end do
          do q = m+1, l-1
            numerator_product = numerator_product * dble(i - q + 1)
          end do
          do q = l+1, eno_order
            numerator_product = numerator_product * dble(i - q + 1)
          end do
          numerator = numerator + numerator_product
        end do

        eno_coeffs(i, j) = eno_coeffs(i, j) + numerator / denominator

      end do

    end do
  end do

end subroutine ENO_Setup

 /*@@
   @routine    ENO_Shutdown
   @date       Mon Apr  8 12:40:44 2002
   @author     Ian Hawke
   @desc 
   Deallocates the coefficient arrays
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine ENO_Shutdown(CCTK_ARGUMENTS)

  USE ENO_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS

  CCTK_INT :: deallocstat

  deallocate(eno_coeffs, STAT = deallocstat)
  if (deallocstat .ne. 0) call CCTK_WARN(0, "Failed to deallocate ENO coefficients.")

end subroutine ENO_Shutdown

 /*@@
   @routine    ENOReconstruct1d
   @date       Sat Apr  6 18:15:29 2002
   @author     Ian Hawke
   @desc 
   Perform a one dimensional reconstruction of a given array using ENO.
   Again, see Shu, section 2.1
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

#define SpaceMask_CheckStateBitsF90_1D(mask,i,type_bits,state_bits) \
  (iand(mask((i)),(type_bits)).eq.(state_bits))

subroutine ENOReconstruct1d(order, nx, v, vminus, vplus, trivial_rp, &
     space_mask, excision_descriptors)
  USE ENO_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: order, nx, i, j, k, r
  CCTK_REAL, dimension(nx) :: v, vplus, vminus
  CCTK_REAL, dimension(order, 1-order:nx+order) :: vdiff

  CCTK_INT, dimension(nx) :: space_mask
  logical, dimension(nx) :: trivial_rp
  logical, dimension(nx) :: excise
  CCTK_INT :: excision_bits, excision_mask
  CCTK_INT, dimension(3) :: excision_descriptors

  CCTK_REAL :: large = 1.d10

  excision_bits=excision_descriptors(1)
  excision_mask=excision_descriptors(2)

  vminus = 0.d0
  vplus = 0.d0
  vdiff = 0.d0

  vdiff(1, 1:nx) = v
  excise = .false.
  trivial_rp = .false.

    
  do i = 1, nx


!!$    Calculate the undivided differences

      do k = 2, order
        do j = max(1, i - order), min(nx, i + order)
          vdiff(k, j) = vdiff(k - 1, j + 1) - vdiff(k - 1, j)
        end do
      end do
      
!!$    Ensure the stencil stays within the grid
      
      vdiff(:, 1 - order : 0) = 1.d10
      vdiff(:, nx + 1 : nx + order) = 1.d10
      
!!$    Find the stencil
      
      r = 0
      
      do j = 2, order
        if ( abs(vdiff(j, i-r-1)) < abs(vdiff(j, i-r)) ) r = r + 1
      end do
      
!!$    Calculate the reconstruction
      
      do j = 0, order - 1
        vminus(i) = vminus(i) + eno_coeffs(r-1, j) * vdiff(1, i-r+j)
        vplus(i) = vplus(i) + eno_coeffs(r, j) * vdiff(1, i-r+j)
      end do

    
  end do

end subroutine ENOReconstruct1d
