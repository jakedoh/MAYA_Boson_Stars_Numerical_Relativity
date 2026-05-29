! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "gauge.h"

subroutine NullEvol_Bdry_whitehole(CCTK_ARGUMENTS)
  implicit none


  CCTK_REAL :: SMass, RR
  CCTK_INT :: i
  !  character(len=500) :: message

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS


  !  call CCTK_INFO(message)
  SMass = .5d0 * null_rwt

  boundary_maskn = 2
  boundary_masks = 2


  do i = 1, 2

     RR = null_rb(i)

     jcn(:,:,i) = (0.,0.)
     jcs(:,:,i) = (0.,0.)

     eth2jcn(:,:,i) = (0.,0.)
     eth2jcs(:,:,i) = (0.,0.)

     ucn(:,:,i) = (0.,0.)
     ucs(:,:,i) = (0.,0.)

#ifdef HORIZON_GAUGE
     bcn(:,:,i) = .5d0 * log(-4.0d0*SMass / cctk_time)
     bcs(:,:,i) = .5d0 * log(-4.0d0*SMass / cctk_time)

     wcn(:,:,i) = (-4.0d0 * SMass / cctk_time *&
          ( 1.0d0 - 2.0d0*SMass / RR) - 1.0d0 ) /RR
     wcs(:,:,i) = (-4.0d0 * SMass / cctk_time *&
          ( 1.0d0 - 2.0d0*SMass / RR) - 1.0d0 ) /RR
#else
     bcn(:,:,i) = 0
     bcs(:,:,i) = 0

     wcn(:,:,i) = -2.0d0 * SMass / RR**2
     wcs(:,:,i) = -2.0d0 * SMass / RR**2
#endif 

     if (first_order_scheme.ne.0) then

        cbcn(:,:,i) = (0.,0.)
        cbcs(:,:,i) = (0.,0.)

        ckcn(:,:,i) = (0.,0.)
        ckcs(:,:,i) = (0.,0.)

        nucn(:,:,i) = (0.,0.)
        nucs(:,:,i) = (0.,0.)

     end if

  end do

end subroutine NullEvol_Bdry_whitehole

subroutine NullEvol_Bdry_flat(CCTK_ARGUMENTS)
  implicit none

  CCTK_INT :: i

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  boundary_maskn = 2 !3
  boundary_masks = 2 !3

  do i = 1, 2 !min(N_radial_pts, 10)

     jcn(:,:,i) = (0.,0.)
     jcs(:,:,i) = (0.,0.)

     eth2jcn(:,:,i) = (0.,0.)
     eth2jcs(:,:,i) = (0.,0.)

     ucn(:,:,i) = (0.,0.)
     ucs(:,:,i) = (0.,0.)

     bcn(:,:,i) = 0
     bcs(:,:,i) = 0

     wcn(:,:,i) = 0
     wcs(:,:,i) = 0

     if (first_order_scheme.ne.0) then

        cbcn(:,:,i) = (0.,0.)
        cbcs(:,:,i) = (0.,0.)

        ckcn(:,:,i) = (0.,0.)
        ckcs(:,:,i) = (0.,0.)

        nucn(:,:,i) = (0.,0.)
        nucs(:,:,i) = (0.,0.)

     end if

  end do

end subroutine NullEvol_Bdry_flat

subroutine NullEvol_Bdry_randomJ(CCTK_ARGUMENTS)
  implicit none

  CCTK_INT :: i

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: random_amp = 1.e-15

  boundary_maskn = 2 !3
  boundary_masks = 2 !3

  ! i = min(N_radial_pts, 10)

  do i = 1, 2 !min(N_radial_pts, 10)

     call random_number(bcn(:,:,:i))
     write(*,*)"random_number:", i, maxval(abs(bcn(:,:,:i))) 
     call random_number(wcn(:,:,:i))
     jcn(:,:,:i) = dcmplx(bcn(:,:,:i), wcn(:,:,:i)) * random_amp

     call random_number(bcs(:,:,:i))
     call random_number(wcs(:,:,:i))
     jcs(:,:,:i) = dcmplx(bcs(:,:,:i), wcs(:,:,:i)) * random_amp

     write(*,*)"RANDOM_J", i, max(maxval(abs(jcn(:,:,:i))), maxval(abs(jcs(:,:,:i)))) 

     eth2jcn(:,:,i) = (0.,0.)
     eth2jcs(:,:,i) = (0.,0.)

     ucn(:,:,i) = (0.,0.)
     ucs(:,:,i) = (0.,0.)

     bcn(:,:,i) = 0
     bcs(:,:,i) = 0

     wcn(:,:,i) = 0
     wcs(:,:,i) = 0

     if (first_order_scheme.ne.0) then

        cbcn(:,:,i) = (0.,0.)
        cbcs(:,:,i) = (0.,0.)

        ckcn(:,:,i) = (0.,0.)
        ckcs(:,:,i) = (0.,0.)

        nucn(:,:,i) = (0.,0.)
        nucs(:,:,i) = (0.,0.)

     end if

  end do

end subroutine NullEvol_Bdry_randomJ
