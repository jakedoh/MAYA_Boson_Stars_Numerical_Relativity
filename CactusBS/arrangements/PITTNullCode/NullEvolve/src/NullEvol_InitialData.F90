! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NullEvol_InitialData(CCTK_ARGUMENTS)
  use NullEvol_sYlm
  use NullEvol_Pulse
  use NullInterp
  implicit none

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)
  
  CCTK_INT :: i, ll, mm, qp, pq
  CCTK_REAL :: PU, qns, pns, i_r
  CCTK_COMPLEX :: jr0, jr1, jr3, jr4, jr5, jr6
  CCTK_REAL :: x_interior, x_exterior, x_mask
  CCTK_INT :: i_interior, i_exterior, i_mask

  CCTK_COMPLEX, dimension(:,:,:) , allocatable :: JCOM
  DECLARE_CCTK_ARGUMENTS
  CCTK_COMPLEX, dimension(null_lsh(1), null_lsh(2),2) :: Ylm
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  call CCTK_INFO("Null Initial Data")

  ! cctk_delta_time set by null step (called below)

  
  eth2jcn = (0.,0.)
  eth2jcs = (0.,0.)

  dxjcn = (0.,0.)
  dxjcs = (0.,0.)

  if (CCTK_EQUALS(initial_J_data,"vanishing_J")) then

     jcn = (0.,0.)
     jcs = (0.,0.)

  else if (CCTK_EQUALS(initial_J_data,"smooth_J")) then

     do i = 1, N_radial_pts
     
        jcn(:,:,i) = (null_xb(i) - 1.d0) / (x_wt(:,:,1) - 1.d0)* j_wt(:,:,1)
        jcs(:,:,i) = (null_xb(i) - 1.d0) / (x_wt(:,:,2) - 1.d0)* j_wt(:,:,2)

     end do

!old J initialization 
!     i_mask = max( maxval(boundary_maskn), maxval(boundary_masks) )
!     do i = 1, N_radial_pts
!        jcn(:,:,i) = ( null_xb(i) - 1.d0 ) / ( null_xb(i_mask) - 1.d0 )*jcn(:,:,i_mask) 
!        jcs(:,:,i) = ( null_xb(i) - 1.d0 ) / ( null_xb(i_mask) - 1.d0 )*jcs(:,:,i_mask) 
!     end do
 
  else if (CCTK_EQUALS(initial_J_data,"vanishing_J_scri")) then

     i_mask = max( maxval(boundary_maskn), maxval(boundary_masks) )

     x_interior = 0.9 * null_xb(i_mask) + 0.1
     x_exterior = 0.1 * null_xb(i_mask) + 0.9

     i_interior = int( ( x_interior - null_xb(1) ) / null_dx ) + 1
     i_exterior = int( ( x_exterior - null_xb(1) ) / null_dx ) + 1

     i_interior = max( i_interior, int( ( null_xb(i_mask) - null_xb(1) ) / null_dx ) + 4 )
     i_exterior = min( i_exterior, N_radial_pts - 3 )

     x_interior = null_xb(i_interior)
     x_exterior = null_xb(i_exterior)

     do i = 1, N_radial_pts

        x_mask =  ( null_xb(i) - x_exterior ) / ( x_interior - x_exterior )
        x_mask = max(0.d0, min(1.d0, x_mask))
        x_mask = (10.d0 + (-15.d0 + 6.d0 * x_mask) * x_mask) * x_mask ** 3

        jcn(:,:,i) = x_mask * jcn(:,:,i)
        jcs(:,:,i) = x_mask * jcs(:,:,i)

     end do

  else if (CCTK_EQUALS(initial_J_data,"rotating_pulse")) then
     call CCTK_INFO("ROTATING PULSE")

     allocate(JCOM(null_lsh(1), null_lsh(2), 2))
     JCOM = 0

     jcn = 0.0
     jcs = 0.0


     do i = 1, N_radial_pts- 1
        if ( null_xb(i) > xmin .and. null_xb(i) < xmax ) then
           do mm = 1, null_lsh(2)
              do ll = 1, null_lsh(1)

                 jcn(ll,mm,i) = WKB_pulse(stereo_q(ll,mm), stereo_p(ll,mm),&
                      null_xb(i), null_rwt, qsmin, qsmax, psmin, psmax,&
                      xmin, xmax, wrot, cctk_time, ID_AMP)

                 if (stereo_q(ll,mm) .ne. 0 .OR. stereo_p(ll,mm) .ne. 0) then
                    qns =  stereo_q(ll,mm) /&
                         (stereo_q(ll,mm)**2 + stereo_p(ll,mm)**2)
                    pns = -stereo_p(ll,mm) /&
                         (stereo_q(ll,mm)**2 + stereo_p(ll,mm)**2)
                 else
                    qns = 1.0d10
                    pns = 1.0d10
                 endif

                 jcs(ll,mm,i) = WKB_pulse(qns, pns,&
                      null_xb(i), null_rwt, qsmin, qsmax, psmin, psmax,&
                      xmin, xmax, wrot, cctk_time, ID_AMP)

              end do
           end do

           call NullInterp_d1(JCOM(:,:,1), jcn(:,:,i), 0_ik, 1_ik)
           call NullInterp_d1(JCOM(:,:,2), jcs(:,:,i), 0_ik, 1_ik)
           call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, JCOM(:,:,1), JCOM(:,:,2), 1_ik)

           call NullInterp_d1(jcn(:,:,i), JCOM(:,:,1), 1_ik, 1_ik)
           call NullInterp_d1(jcs(:,:,i), JCOM(:,:,2), 1_ik, 1_ik)
           call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, jcn(:,:,i), jcs(:,:,i), 2_ik)

        endif

     end do
     deallocate(JCOM)

  else if (CCTK_EQUALS(initial_J_data,"compact_J")) then

     if (use_rsYlm.ne.0) then
        call CCTK_INFO("Using real spherical Harmonics")
        call rsYlm(2_ik,ID_l, ID_m, null_lsh(1), null_lsh(2), zeta, Ylm)
     else
        call CCTK_INFO("Using standard spherical Harmonics")
        call sYlm(2_ik,ID_l, ID_m, null_lsh(1), null_lsh(2), zeta, Ylm)
     endif

     do i = 1, N_radial_pts
        if ( null_xb(i) > xmin .and. null_xb(i) < xmax ) then
           PU = 2.0d0**(2*ID_power) * &
                ((null_xb(i) - xmin) * (xmax - null_xb(i)))**ID_power &
                / ( xmax - xmin )**(2*ID_power)
        else 
           PU = 0
        end if

        ! remember jcn(:,:,i) is a larger array than zeta(:,:)

        jcn(1:null_lsh(1),1:null_lsh(2),i) =&
             ID_AMP * PU *Ylm(:,:,1)
        jcs(1:null_lsh(1),1:null_lsh(2),i) =&
             ID_AMP * PU *Ylm(:,:,2)

        if (ID_l .eq. 2 .AND. ID_m .eq. 0) then
           jcn(1:null_lsh(1),1:null_lsh(2),i) =&
                ID_AMP * PU * 4.0*(zeta*zeta)/(1+zeta*conjg(zeta))**2
           jcs(1:null_lsh(1),1:null_lsh(2),i) =&
                ID_AMP * PU * 4.0*(zeta*zeta)/(1+zeta*conjg(zeta))**2
        endif

     end do

  else if (CCTK_EQUALS(initial_J_data,"extracted_J")) then
     ! already done in extraction

  else if (CCTK_EQUALS(initial_J_data,"polynomial_J")) then

     do i = 1, N_radial_pts
!J=(x-1)^3/(x_E-1)^3[J_wt+(x-x_E)(J,x_wt-3J_wt/(x_E)-1)]
!J_x_wt = exp(2beta_wt)rwt/(1-x_wt)^2
                jcn(:,:,i) = (null_xb(i) - 1.d0)**3 / (x_wt(:,:,1) - 1)**3&
                           * (j_wt(:,:,1) + (null_xb(i) - x_wt(:,:,1))&
                           * (j_l(:,:,1) * exp(2.d0 * beta_wt(:,:,1))&
                           * null_rwt / (1 - x_wt(:,:,1))**2&
                           - 3.d0 * j_wt(:,:,1) / (x_wt(:,:,1) - 1.d0)))
                jcs(:,:,i) = (null_xb(i) - 1.d0)**3 / (x_wt(:,:,2) - 1)**3&
                           * (j_wt(:,:,2) + (null_xb(i) - x_wt(:,:,2))&
                           * (j_l(:,:,2) * exp(2.d0 * beta_wt(:,:,2))&
                           * null_rwt / (1 - x_wt(:,:,2))**2&
                           - 3.d0 * j_wt(:,:,2) / (x_wt(:,:,2) - 1.d0)))

     end do
 
  else if (CCTK_EQUALS(initial_J_data,"fitted_linearized_J")) then
     ! use a special solution adapted to a previous run
    
     ! FORMULA FOR SETTING (NON-ZERO) INITIAL DATA
     ! J=j2(r) Y_{2,2} + j2(r)^* Y_{2,-2} + j0(r) Y_{2,0}
     ! where
     ! j2(r) = (Jcoeff_r0r+i*Jcoeff_r0i) + (Jcoeff_r1r+i*Jcoeff_r1i)/r
     !         +(Jcoeff_r3r+i*Jcoeff_r3i)/r^3
     ! j0(r)= Jcoeff_r0 + Jcoeff_r1/r + Jcoeff_r3/r^3
     !    with ^* meaning complex conjugate. 
     
     ! i_r = (1-x)/(x*rwt)
     
     call CCTK_INFO("Using fitted linearized initial data for J")
     
     jr0 = dcmplx(Jcoeff_r0r, Jcoeff_r0i)
     jr1 = dcmplx(Jcoeff_r1r, Jcoeff_r1i)
     jr3 = dcmplx(Jcoeff_r3r, Jcoeff_r3i)
     jr4 = Jcoeff_r0
     jr5 = Jcoeff_r1
     jr6 = Jcoeff_r3

     ! Y22 term
     call sYlm(2_ik,2_ik,2_ik, null_lsh(1), null_lsh(2), zeta, Ylm)
     do i = 1, N_radial_pts
        i_r = (1.0-null_xb(i))/(null_xb(i)*null_rwt)
        jcn(:,:,i) = (jr0 + jr1*i_r + jr3*i_r**3) * Ylm(:,:,1)
        jcs(:,:,i) = (jr0 + jr1*i_r + jr3*i_r**3) * Ylm(:,:,2)
     end do
     ! Y2-2 term
     call sYlm(2_ik,2_ik,-2_ik, null_lsh(1), null_lsh(2), zeta, Ylm)
     do i = 1, N_radial_pts
        i_r = (1.0-null_xb(i))/(null_xb(i)*null_rwt)
        jcn(:,:,i) = jcn(:,:,i) + conjg((jr0 + jr1*i_r + jr3*i_r**3)) * Ylm(:,:,1)
        jcs(:,:,i) = jcs(:,:,i) + conjg((jr0 + jr1*i_r + jr3*i_r**3)) * Ylm(:,:,2)
     end do
     ! Y20 term
     call sYlm(2_ik,2_ik,0_ik, null_lsh(1), null_lsh(2), zeta, Ylm)
     do i = 1, N_radial_pts
        i_r = (1.0-null_xb(i))/(null_xb(i)*null_rwt)
        jcn(:,:,i) = jcn(:,:,i) + (jr4 + jr5*i_r + jr6*i_r**3) * Ylm(:,:,1)
        jcs(:,:,i) = jcs(:,:,i) + (jr4 + jr5*i_r + jr6*i_r**3) * Ylm(:,:,2)
     end do
     
  end if

!mask J and eth2_J initially

  do i = 1, N_radial_pts

      jcn(:,:,i) = EG_mask * jcn(:,:,i)
      jcs(:,:,i) = EG_mask * jcs(:,:,i)

      eth2jcn(:,:,i) = EQ_mask * eth2jcn(:,:,i)
      eth2jcs(:,:,i) = EQ_mask * eth2jcs(:,:,i)

  end do

!initialize dx_J

  dxjcn(:,:,1) = -0.5 * (3.*jcn(:,:,1) - 4.*jcn(:,:,2) + jcn(:,:,3)) / null_dx
  dxjcs(:,:,1) = -0.5 * (3.*jcs(:,:,1) - 4.*jcs(:,:,2) + jcs(:,:,3)) / null_dx

  do i = 2, N_radial_pts-1

     dxjcn(:,:,i) = 0.5 * (jcn(:,:,i+1) - jcn(:,:,i-1)) / null_dx
     dxjcs(:,:,i) = 0.5 * (jcs(:,:,i+1) - jcs(:,:,i-1)) / null_dx

  end do

  dxjcn(:,:,N_radial_pts) = 0.5 * (3.*jcn(:,:,N_radial_pts) &
              - 4.*jcn(:,:,N_radial_pts-1) + jcn(:,:,N_radial_pts-2)) / null_dx
  dxjcs(:,:,N_radial_pts) = 0.5 * (3.*jcs(:,:,N_radial_pts) &
              - 4.*jcs(:,:,N_radial_pts-1) + jcs(:,:,N_radial_pts-2)) / null_dx

!initialize array for the radial profile
  qp = int((null_lsh(1)-1)/2) + 1
  pq = int((null_lsh(2)-1)/2) + 1

  jcn_rad(:) = jcn(qp,pq,:)
  jcs_rad(:) = jcs(qp,pq,:)

  dxjcn_rad(:) = dxjcn(qp,pq,:)
  dxjcs_rad(:) = dxjcs(qp,pq,:)

!  write (*,*) 'NULL_INITIALDATA: J_scri N are: ', maxval(abs(jcn(:,:,N_radial_pts)))
!  write (*,*) 'NULL_INITIALDATA: J_scri S are: ', maxval(abs(jcs(:,:,N_radial_pts)))

end subroutine NullEvol_InitialData


subroutine NullEvol_InitializeArrays(CCTK_ARGUMENTS)
  implicit none


  CCTK_COMPLEX :: ii = (0., 1.)

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  jcs = ii*1.d+50; jcs_p = ii*1.d+50; jcn = ii*1.d+50; jcn_p = ii*1.d+50
  eth2jcs = ii*1.d+50; eth2jcs_p = ii*1.d+50; eth2jcn = ii*1.d+50; eth2jcn_p = ii*1.d+50
  ucs = ii*1.d+50; ucs_p = ii*1.d+50; ucn = ii*1.d+50; ucn_p = ii*1.d+50
  wcs = 1.d+50; wcs_p = 1.d+50; wcn = 1.d+50; wcn_p = 1.d+50
  bcs = 1.d+50; bcs_p = 1.d+50; bcn = 1.d+50; bcn_p = 1.d+50

  if (first_order_scheme.ne.0) then

     nucs = ii*1.d+50; nucs_p = ii*1.d+50; nucn = ii*1.d+50; nucn_p = ii*1.d+50
     ckcs = ii*1.d+50; ckcs_p = ii*1.d+50; ckcn = ii*1.d+50; ckcn_p = ii*1.d+50
     cbcs = ii*1.d+50; cbcs_p = ii*1.d+50; cbcn = ii*1.d+50; cbcn_p = ii*1.d+50

  end if
  
  jcn_rad = ii*1.d+50; jcs_rad = ii*1.d+50
  dxjcn_rad = ii*1.d+50; dxjcs_rad = ii*1.d+50

end subroutine NullEvol_InitializeArrays

