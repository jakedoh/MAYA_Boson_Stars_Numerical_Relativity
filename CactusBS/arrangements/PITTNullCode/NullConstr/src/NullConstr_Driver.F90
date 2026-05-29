! vim: syntax=fortran

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#define nq null_lsh(1)
#define np null_lsh(2)

#define stereo .true.

subroutine NullConstr_Driver(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars
  use NullInterp
  use NullEvol_Mask
  use NullConstr_Util
  use NullConstr_R00
  use NullConstr_R01
  use NullConstr_R0A

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i

  CCTK_REAL, allocatable, save, dimension(:) :: &
       dr_dx,d2r_dx2,dr_dxh,d2r_dx2h

  CCTK_REAL, allocatable, save, dimension(:,:) :: &
       w_00,w_01,w_04,w_11,w_14,&
       b_00,b_01,b_04,b_11,b_14,&
       k_00,k_01,k_04,k_11,k_14

  CCTK_REAL, allocatable, save, dimension(:,:,:) :: kc, kc_p

  CCTK_COMPLEX, allocatable, save, dimension(:,:) :: &
       w_02,w_03,w_12,w_13,w_22,w_23,w_24,w_33,w_34,&
       b_02,b_03,b_12,b_13,b_22,b_23,b_24,b_33,b_34,&
       k_02,k_03,k_12,k_13,k_22,k_23,k_24,k_33,k_34,&
       j_00,j_01,j_04,j_11,j_14,&
       j_02,j_03,j_12,j_13,j_22,j_23,j_24,j_33,j_34,&
       jb_00,jb_01,jb_04,jb_11,jb_14,&
       jb_02,jb_03,jb_12,jb_13,jb_22,jb_23,jb_24,jb_33,jb_34,&
       u_00,u_01,u_04,u_11,u_14,&
       u_02,u_03,u_12,u_13,u_22,u_23,u_24,u_33,u_34,&
       ub_00,ub_01,ub_04,ub_11,ub_14,&
       ub_02,ub_03,ub_12,ub_13,ub_22,ub_23,ub_24,ub_33,ub_34

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  LOGICAL, save :: FirstTime = .true.
  integer :: l1, l2, reduce_handle
  CCTK_INT :: mn, reval, global_mn

  call CCTK_ReductionArrayHandle(reduce_handle, "minimum");
  if (reduce_handle .lt. 0 ) then
     call CCTK_WARN(0,"Could not get reduction handle")
  endif

  if (stereo) then
     mn = min(minval(boundary_masks(1:nq,1:np)),&
          minval(boundary_maskn(1:nq,1:np)))
  else
     mn = minval(boundary_maskn(1:nq,1:np))
  end if

  call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
       mn, global_mn, CCTK_VARIABLE_INT)

  if (reval .ne. 0 ) then
     call CCTK_WARN(0,"Error in obtaining Minimum of mask")
  endif

  call CCTK_ActiveTimeLevels(l1, cctkGH, "NullVars::realcharfuncs")
  call CCTK_ActiveTimeLevels(l2, cctkGH, "NullVars::cmplxcharfuncs_basic")

  if (min(l1,l2).lt.2) then
     call CCTK_WARN(1, "cannot calculate constraints -- not enough allocated time-levels")
     !     Null_R00 = 1.e+10
     return
  end if
  if (FirstTime) then
     FirstTime = .false.
     allocate(dr_dx(N_radial_pts),d2r_dx2(N_radial_pts),dr_dxh(N_radial_pts),d2r_dx2h(N_radial_pts),&
          kc(nq,np,N_radial_pts),kc_p(nq,np,N_radial_pts), &
          w_00(nq,np),w_01(nq,np),&
          w_04(nq,np),w_11(nq,np),&
          w_14(nq,np), w_02(nq,np),&
          w_03(nq,np),w_12(nq,np),&
          w_13(nq,np),w_22(nq,np),&
          w_23(nq,np),w_24(nq,np),&
          w_33(nq,np),w_34(nq,np),&
          b_00(nq,np),b_01(nq,np),&
          b_04(nq,np),b_11(nq,np),&
          b_14(nq,np),b_02(nq,np),&
          b_03(nq,np),b_12(nq,np),&
          b_13(nq,np),b_22(nq,np),&
          b_23(nq,np),b_24(nq,np),&
          b_33(nq,np),b_34(nq,np),&
          k_00(nq,np),k_01(nq,np),&
          k_04(nq,np),k_11(nq,np),&
          k_14(nq,np),k_02(nq,np),&
          k_03(nq,np),k_12(nq,np),&
          k_13(nq,np),k_22(nq,np),&
          k_23(nq,np),k_24(nq,np),&
          k_33(nq,np),k_34(nq,np),&
          j_00(nq,np),j_01(nq,np),&
          j_04(nq,np),j_11(nq,np),&
          j_14(nq,np),j_02(nq,np),&
          j_03(nq,np),j_12(nq,np),&
          j_13(nq,np),j_22(nq,np),&
          j_23(nq,np),j_24(nq,np),&
          j_33(nq,np),j_34(nq,np),&
          jb_00(nq,np),jb_01(nq,np),&
          jb_04(nq,np),jb_11(nq,np),&
          jb_14(nq,np),jb_02(nq,np),&
          jb_03(nq,np),jb_12(nq,np),&
          jb_13(nq,np),jb_22(nq,np),&
          jb_23(nq,np),jb_24(nq,np),&
          jb_33(nq,np),jb_34(nq,np),&
          u_00(nq,np),u_01(nq,np),&
          u_04(nq,np),u_11(nq,np),&
          u_14(nq,np),u_02(nq,np),&
          u_03(nq,np),u_12(nq,np),&
          u_13(nq,np),u_22(nq,np),&
          u_23(nq,np),u_24(nq,np),&
          u_33(nq,np),u_34(nq,np),&
          ub_00(nq,np),ub_01(nq,np),&
          ub_04(nq,np),ub_11(nq,np),&
          ub_14(nq,np),ub_02(nq,np),&
          ub_03(nq,np),ub_12(nq,np),&
          ub_13(nq,np),ub_22(nq,np),&
          ub_23(nq,np),ub_24(nq,np),&
          ub_33(nq,np),ub_34(nq,np))
  end if

  Null_R00 = 0; Null_R01 = 0; Null_R0A = 0
  if (cctk_iteration .eq. 0) then
    ! Can't calculate time derivatives at t = 0
    return
  endif

  call x2r_derivs(dr_dx,d2r_dx2,dr_dxh,d2r_dx2h,N_radial_pts,xb,xbh)

  kc = dsqrt(1+dble(jcn*conjg(jcn)))
  kc_p = dsqrt(1+dble(jcn_p*conjg(jcn_p)))


 ! North patch:

  do i = global_mn+1, N_radial_pts  

     call real_derivs (wcn,wcn_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,dr_dx(i),d2r_dx2(i), &
          w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
          w_22,w_23,w_24,w_33,w_34)
     call real_derivs (bcn,bcn_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,dr_dx(i),d2r_dx2(i), &
          b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
          b_22,b_23,b_24,b_33,b_34)
     call real_derivs (kc,kc_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,dr_dx(i),d2r_dx2(i), &
          k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
          k_22,k_23,k_24,k_33,k_34)
     call cmplx_derivs &
          (jcn,jcn_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,2_ik,dr_dx(i),d2r_dx2(i),.FALSE., &
          j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
          j_22,j_23,j_24,j_33,j_34, &
          jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
          jb_22,jb_23,jb_24,jb_33,jb_34)
     call cmplx_derivs &
          (ucn,ucn_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,1_ik,dr_dx(i),d2r_dx2(i),.TRUE., &
          u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
          u_22,u_23,u_24,u_33,u_34, &
          ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
          ub_22,ub_23,ub_24,ub_33,ub_34)

     !calculate R00n,R01n,R0An

     call NullConstr_R00_calc (Null_R00(:,:,i),nq,np,null_rb(i), &
          w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
          w_22,w_23,w_24,w_33,w_34, &
          b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
          b_22,b_23,b_24,b_33,b_34, &
          k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
          k_22,k_23,k_24,k_33,k_34, &
          j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
          j_22,j_23,j_24,j_33,j_34, &
          jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
          jb_22,jb_23,jb_24,jb_33,jb_34, &
          u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
          u_22,u_23,u_24,u_33,u_34, &
          ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
          ub_22,ub_23,ub_24,ub_33,ub_34)

     call NullConstr_R0A_calc (Null_R0A(:,:,i),nq,np,null_rb(i), &
          w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
          w_22,w_23,w_24,w_33,w_34, &
          b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
          b_22,b_23,b_24,b_33,b_34, &
          k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
          k_22,k_23,k_24,k_33,k_34, &
          j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
          j_22,j_23,j_24,j_33,j_34, &
          jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
          jb_22,jb_23,jb_24,jb_33,jb_34, &
          u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
          u_22,u_23,u_24,u_33,u_34, &
          ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
          ub_22,ub_23,ub_24,ub_33,ub_34)

     call NullConstr_R01_calc (Null_R01(:,:,i),nq,np,null_rb(i), &
          w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
          w_22,w_23,w_24,w_33,w_34, &
          b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
          b_22,b_23,b_24,b_33,b_34, &
          k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
          k_22,k_23,k_24,k_33,k_34, &
          j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
          j_22,j_23,j_24,j_33,j_34, &
          jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
          jb_22,jb_23,jb_24,jb_33,jb_34, &
          u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
          u_22,u_23,u_24,u_33,u_34, &
          ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
          ub_22,ub_23,ub_24,ub_33,ub_34)
!         ub_22,ub_23,ub_24,ub_33,ub_34,omm,cctk_time,Ylm_0,1)

     call NullEvol_remask (i, boundary_maskn, evolution_maskn)

     Null_R00(:,:,i) = evolution_maskn * Null_R00(:,:,i)
     Null_R0A(:,:,i) = evolution_maskn * Null_R0A(:,:,i)
     Null_R01(:,:,i) = evolution_maskn * Null_R01(:,:,i)

 end do
 Null_R0A_r = dble(Null_R0A)
 Null_R0A_i = dimag(Null_R0A)

 ! south patch (only if stereographic)

 if (stereo) then

    kc = dsqrt(1+dble(jcs*conjg(jcs)))
    kc_p = dsqrt(1+dble(jcs_p*conjg(jcs_p)))

    Null_R00_south = 0; Null_R01_south = 0; Null_R0A_south = 0

    do i = global_mn+1, N_radial_pts  

        call real_derivs (wcs,wcs_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,dr_dx(i),d2r_dx2(i), &
             w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
             w_22,w_23,w_24,w_33,w_34)
        call real_derivs (bcs,bcs_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,dr_dx(i),d2r_dx2(i), &
             b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
             b_22,b_23,b_24,b_33,b_34)
        call real_derivs (kc,kc_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,dr_dx(i),d2r_dx2(i), &
             k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
             k_22,k_23,k_24,k_33,k_34)
        call cmplx_derivs &
             (jcs,jcs_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,2_ik,dr_dx(i),d2r_dx2(i),.FALSE., &
             j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
             j_22,j_23,j_24,j_33,j_34, &
             jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
             jb_22,jb_23,jb_24,jb_33,jb_34)
        call cmplx_derivs &
             (ucs,ucs_p,cctk_delta_time,dx,N_radial_pts,null_lsh,i,1_ik,dr_dx(i),d2r_dx2(i),.TRUE., &
             u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
             u_22,u_23,u_24,u_33,u_34, &
             ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
             ub_22,ub_23,ub_24,ub_33,ub_34)

        !calculate R00s,R01s,R0As
        call NullConstr_R00_calc (Null_R00_south(:,:,i),nq,np,null_rb(i), &
             w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
             w_22,w_23,w_24,w_33,w_34, &
             b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
             b_22,b_23,b_24,b_33,b_34, &
             k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
             k_22,k_23,k_24,k_33,k_34, &
             j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
             j_22,j_23,j_24,j_33,j_34, &
             jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
             jb_22,jb_23,jb_24,jb_33,jb_34, &
             u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
             u_22,u_23,u_24,u_33,u_34, &
             ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
             ub_22,ub_23,ub_24,ub_33,ub_34)

        call NullConstr_R0A_calc (Null_R0A_south(:,:,i),nq,np,null_rb(i), &
             w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
             w_22,w_23,w_24,w_33,w_34, &
             b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
             b_22,b_23,b_24,b_33,b_34, &
             k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
             k_22,k_23,k_24,k_33,k_34, &
             j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
             j_22,j_23,j_24,j_33,j_34, &
             jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
             jb_22,jb_23,jb_24,jb_33,jb_34, &
             u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
             u_22,u_23,u_24,u_33,u_34, &
             ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
             ub_22,ub_23,ub_24,ub_33,ub_34)

        call NullConstr_R01_calc (Null_R01_south(:,:,i),nq,np,null_rb(i), &
             w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
             w_22,w_23,w_24,w_33,w_34, &
             b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
             b_22,b_23,b_24,b_33,b_34, &
             k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
             k_22,k_23,k_24,k_33,k_34, &
             j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
             j_22,j_23,j_24,j_33,j_34, &
             jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
             jb_22,jb_23,jb_24,jb_33,jb_34, &
             u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
             u_22,u_23,u_24,u_33,u_34, &
             ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
             ub_22,ub_23,ub_24,ub_33,ub_34)
!            ub_22,ub_23,ub_24,ub_33,ub_34,omm,cctk_time,Ylm_0,2)

        call NullEvol_remask (i, boundary_masks, evolution_masks)

        Null_R00_south(:,:,i) = evolution_masks * Null_R00_south(:,:,i)
        Null_R0A_south(:,:,i) = evolution_masks * Null_R0A_south(:,:,i)
        Null_R01_south(:,:,i) = evolution_masks * Null_R01_south(:,:,i)

        call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, Null_R0A(:,:,i),Null_R0A_south(:,:,i), 1_ik)
        call NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, Null_R00(:,:,i), Null_R00_south(:,:,i))
        call NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, Null_R01(:,:,i), Null_R01_south(:,:,i))

     end do

  end if

  call CCTK_INFO("CALCULATED NULL_R00, NULL_R0A, NULL_R01")
end subroutine NullConstr_Driver
