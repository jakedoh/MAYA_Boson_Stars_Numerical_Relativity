! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NoExcision_CGInit_1 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS_NoExcision_CGInit_1
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: my_level, n_levels

  nx = cctk_lsh(1); ny = cctk_lsh(2); nz = cctk_lsh(3)

  loop_control = 0

  call NoExcision_levelinfo ( cctkGH, my_level, n_levels )

  if ( my_level == n_levels-1 ) then

    where (abs(gxx) < 1d-8 .and. abs(alp) < 1d-8)
      nes_mask = 1
    elsewhere
      nes_mask = 0
    end where

    call CCTK_INFO ( 'Starting smoothing procedure' )

! This set of OpenMP directives have been removed due to triggering an
! internal compiler error in the Intel 15.0.0 compiler. As this section
! of the code is only compiled once, it shouldn't affect performance in
! any significant way even for compilers that are able to generate
! parallel code (none of the Intel compilers do this anyway).

!!$OMP PARALLEL WORKSHARE

    resgxx = zero; resgxy = zero; resgxz = zero
    resgyy = zero; resgyz = zero; resgzz = zero
    reskxx = zero; reskxy = zero; reskxz = zero
    reskyy = zero; reskyz = zero; reskzz = zero
    res = zero; resx = zero; resy = zero; resz = zero

    dgxx = zero; dgxy = zero; dgxz = zero
    dgyy = zero; dgyz = zero; dgzz = zero
    dkxx = zero; dkxy = zero; dkxz = zero
    dkyy = zero; dkyz = zero; dkzz = zero
    d = zero; dx = zero; dy = zero; dz = zero

    qgxx = zero; qgxy = zero; qgxz = zero
    qgyy = zero; qgyz = zero; qgzz = zero
    qkxx = zero; qkxy = zero; qkxz = zero
    qkyy = zero; qkyz = zero; qkzz = zero
    q = zero; qx = zero; qy = zero; qz = zero

    redgxx = zero; redgxy = zero; redgxz = zero
    redgyy = zero; redgyz = zero; redgzz = zero
    redkxx = zero; redkxy = zero; redkxz = zero
    redkyy = zero; redkyz = zero; redkzz = zero
    red = zero; redx = zero; redy = zero; redz = zero

!!$OMP END PARALLEL WORKSHARE

    ! r = b - A x.
    ! Since x=0 and we actually use A':   b = -A' 0 and  r = b = -A' 0.

    call residual_all ( gxx, gxy, gxz, gyy, gyz, gzz, &
                        kxx, kxy, kxz, kyy, kyz, kzz, &
                        alp, betax, betay, betaz, &
                        resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                        reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                        res, resx, resy, resz, &
                        nes_mask, -one, smoothing_order )

    ! d = r = b.

    call residual_all ( gxx, gxy, gxz, gyy, gyz, gzz, &
                        kxx, kxy, kxz, kyy, kyz, kzz, &
                        alp, betax, betay, betaz, &
                        dgxx, dgxy, dgxz, dgyy, dgyz, dgzz, &
                        dkxx, dkxy, dkxz, dkyy, dkyz, dkzz, &
                        d, dx, dy, dz, &
                        nes_mask, -one, smoothing_order )

    ! red = r*r.

    call multiply ( resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                    reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                    res, resx, resy, resz, &
                    resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                    reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                    res, resx, resy, resz, &
                    redgxx, redgxy, redgxz, redgyy, redgyz, redgzz, &
                    redkxx, redkxy, redkxz, redkyy, redkyz, redkzz, &
                    red, redx, redy, redz, nes_mask, red_mask, &
                    .true., smoothing_order )
                     
    call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
    if ( sum_handle .lt. 0 ) then
      call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
    end if

    call CCTK_ReductionArrayHandle ( infnorm_handle, 'norm_inf' )
    if ( infnorm_handle .lt. 0 ) then
      call CCTK_WARN(0,'Could not obtain a handle for norm_inf reduction')
    end if

    loop_counter = 0
    loop_control = 1

    sym_selector = 1

  end if

end subroutine NoExcision_CGInit_1


subroutine NoExcision_CGInit_2 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: ierr
  integer :: i
  character(len=56) :: conv_message

  if ( loop_control == 1 ) then
    ! delta_new = r^T r.
    call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, lsumred, &
                                        delta_new, 16, CCTK_VARIABLE_REAL)
    if ( ierr < 0 ) call CCTK_WARN ( 0, 'Could not perform reduction of local 1D array' )

    call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, infnorm_handle, &
                                        linfred, infnormresid, 16, &
                                        CCTK_VARIABLE_REAL)
    if ( ierr < 0 ) call CCTK_WARN ( 0, 'Could not perform reduction of local 1D array' )

    where ( cont ) infnormresid = sqrt(infnormresid)

    ! Check if some variables have already converged. This happens when the
    ! variable is identically zero.

    do i = 1, 16 
      if ( cont(i) ) then
        if ( infnormresid(i) < smoothing_eps ) then
          write ( conv_message, '(a23,i8,a20,a5)' ) 'CG method converged in ', &
                loop_counter, ' steps for variable ', var_names(i)
          call CCTK_INFO ( conv_message )
          cont(i) = .false.
        end if 
      end if
    end do

    if ( .not. any ( cont ) ) then
      loop_control = 0
    end if

  end if

end subroutine NoExcision_CGInit_2
    

subroutine NoExcision_CG_1 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  nx = cctk_lsh(1); ny=cctk_lsh(2); nz = cctk_lsh(3)

  loop_counter = loop_counter + 1

  ! Since d is zero outside of the active region we have:
  ! q = A d = A' d.

  call residual_all ( dgxx, dgxy, dgxz, dgyy, dgyz, dgzz, &
                      dkxx, dkxy, dkxz, dkyy, dkyz, dkzz, &
                      d, dx, dy, dz, &
                      qgxx, qgxy, qgxz, qgyy, qgyz, qgzz, &
                      qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, &
                      q, qx, qy, qz, &
                      nes_mask, one, smoothing_order )

  ! red = d*q = d*(A d)

  call multiply ( dgxx, dgxy, dgxz, dgyy, dgyz, dgzz, &
                  dkxx, dkxy, dkxz, dkyy, dkyz, dkzz, &
                  d, dx, dy, dz, &
                  qgxx, qgxy, qgxz, qgyy, qgyz, qgzz, &
                  qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, &
                  q, qx, qy, qz, &
                  redgxx, redgxy, redgxz, redgyy, redgyz, redgzz, &
                  redkxx, redkxy, redkxz, redkyy, redkyz, redkzz, &
                  red, redx, redy, redz, nes_mask, red_mask,&
                  .false., smoothing_order )

  sym_selector = 2

end subroutine NoExcision_CG_1


subroutine NoExcision_CG_2 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  if ( loop_control == 1 ) then
    call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, lsumred, &
                                        alpha, 16, CCTK_VARIABLE_REAL)
    if ( ierr < 0 ) call CCTK_WARN ( 0, 'Could not perform reduction of local 1D array' )

    ! alpha = delta_new / ( d^T A d ).
    where ( cont ) alpha = delta_new / alpha
  end if

end subroutine NoExcision_CG_2


subroutine NoExcision_CG_3 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer :: i

  nx = cctk_lsh(1); ny=cctk_lsh(2); nz = cctk_lsh(3)

  ! x = x + alpha * d.

  call multiply_sum ( gxx, gxy, gxz, gyy, gyz, gzz, &
                      kxx, kxy, kxz, kyy, kyz, kzz, &
                      alp, betax, betay, betaz, &
                      dgxx, dgxy, dgxz, dgyy, dgyz, dgzz, &
                      dkxx, dkxy, dkxz, dkyy, dkyz, dkzz, &
                      d, dx, dy, dz, &
                      (/ (one, i=1, 16) /), alpha, nes_mask )
                      
  if ( mod ( loop_counter, 50 ) == 0 ) then

    ! Restart:
    ! r = b - A x = - A' x

    call residual_all ( gxx, gxy, gxz, gyy, gyz, gzz, &
                        kxx, kxy, kxz, kyy, kyz, kzz, &
                        alp, betax, betay, betaz, &
                        resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                        reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                        res, resx, resy, resz, &
                        nes_mask, -one, smoothing_order )
  else

    ! r = r - alpha q.

    call multiply_sum ( resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                        reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                        res, resx, resy, resz, &
                        qgxx, qgxy, qgxz, qgyy, qgyz, qgzz, &
                        qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, &
                        q, qx, qy, qz, &
                        (/ (one, i=1, 16) /), -alpha, nes_mask )
                      
  end if

  ! delta_old = delta_new.

  where ( cont ) delta_old = delta_new

  ! red = r*r.
  call multiply ( resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                  reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                  res, resx, resy, resz, &
                  resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                  reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                  res, resx, resy, resz, &
                  redgxx, redgxy, redgxz, redgyy, redgyz, redgzz, &
                  redkxx, redkxy, redkxz, redkyy, redkyz, redkzz, &
                  red, redx, redy, redz, nes_mask, red_mask, &
                  .true., smoothing_order )

  sym_selector = 3

end subroutine NoExcision_CG_3


subroutine NoExcision_CG_4 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  if ( loop_control == 1 ) then
    ! delta_new = r^T r.
    call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, lsumred, &
                                        delta_new, 16, CCTK_VARIABLE_REAL)
    if ( ierr < 0 ) call CCTK_WARN ( 0, 'Could not perform reduction of local 1D array' )

    call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, infnorm_handle, &
                                        linfred, infnormresid, 16, &
                                        CCTK_VARIABLE_REAL)
    if ( ierr < 0 ) call CCTK_WARN ( 0, 'Could not perform reduction of local 1D array' )

    where ( cont ) infnormresid = sqrt(infnormresid)
  end if


end subroutine NoExcision_CG_4


subroutine NoExcision_CG_5 (CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer :: i
  character(len=192) :: res_message, var_message
  character(len=56) :: conv_message
  character(len=18) :: iter_message

  nx = cctk_lsh(1); ny=cctk_lsh(2); nz = cctk_lsh(3)

  ! beta = delta_new / delta_old

  where ( cont ) beta = delta_new / delta_old

  ! d = r + beta d

  call multiply_sum ( dgxx, dgxy, dgxz, dgyy, dgyz, dgzz, &
                      dkxx, dkxy, dkxz, dkyy, dkyz, dkzz, &
                      d, dx, dy, dz, &
                      resgxx, resgxy, resgxz, resgyy, resgyz, resgzz, &
                      reskxx, reskxy, reskxz, reskyy, reskyz, reskzz, &
                      res, resx, resy, resz, &
                      beta, (/ (one, i=1, 16) /), nes_mask )
  
  if ( verbose > 0 .and. mod(loop_counter,100) == 0 ) then
    write (var_message, '(16(a3,a9))' ) (' | ', var_names(i), i=1, 16)
    call CCTK_INFO ( var_message )
  end if
  if ( verbose > 0 .and. mod(loop_counter,10) == 0 ) then
    write (iter_message, '(a10,i8)' ) 'Iteration ', loop_counter
    write (res_message, '(16(a3,es9.3))' ) (' | ', infnormresid(i), i=1, 16 )
    call CCTK_INFO ( iter_message )
    call CCTK_INFO ( res_message )
  end if
  
  ! Check if any variables have converged during this iteration.

  do i = 1, 16
    if ( cont(i) ) then
      if ( infnormresid(i) < smoothing_eps ) then
        write ( conv_message, '(a23,i8,a20,a5)' ) 'CG method converged in ', &
                loop_counter, ' steps for variable ', var_names(i)
        call CCTK_INFO ( conv_message )
        cont(i) = .false.
      end if
    end if
  end do

  ! If all variables have converged we exit.

  if ( .not. any ( cont ) ) then
    loop_control = 0
  end if

  sym_selector = 4

end subroutine NoExcision_CG_5


! I leave this in here in case somebody uses CartGrid3D symmetries...
! Not tested, though...

subroutine NoExcision_SetSym(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS

  CCTK_INT :: ierr
  CCTK_INT, dimension(3) :: sym

  sym = 1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resgxx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resgyy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resgzz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::reskxx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::reskyy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::reskzz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dgxx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dgyy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dgzz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dkxx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dkyy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dkzz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qgxx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qgyy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qgzz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qkxx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qkyy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qkzz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::res' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::d' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::q' )

  call SetCartSymGN ( ierr, cctkGH, sym, 'noexcision::cg_red_all' )

  sym(1) = -1; sym(2) = -1; sym(3) = 1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resgxy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::reskxy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dgxy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dkxy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qgxy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qkxy' )

  sym(1) = -1; sym(2) = 1; sym(3) = -1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resgxz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::reskxz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dgxz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dkxz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qgxz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qkxz' )

  sym(1) = 1; sym(2) = -1; sym(3) = -1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resgyz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::reskyz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dgyz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dkyz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qgyz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qkyz' )

  sym(1) = -1; sym(2) = 1; sym(3) = 1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dx' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qx' )

  sym(1) = 1; sym(2) = -1; sym(3) = 1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dy' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qy' )

  sym(1) = 1; sym(2) = 1; sym(3) = -1

  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::resz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::dz' )
  call SetCartSymVN ( ierr, cctkGH, sym, 'noexcision::qz' )

end subroutine NoExcision_SetSym


subroutine NoExcision_CGApplySym(CCTK_ARGUMENTS)

  use NoExcision_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer,  parameter :: ik = kind (izero)
  CCTK_INT :: ierr

  if ( loop_control > 0 ) then

    if ( sym_selector <1 .or. sym_selector > 4 ) then
      call CCTK_WARN ( 0, 'Internal error. Inconsistent symmetry selector' )
    end if
  
    select case ( sym_selector )
    case (1, 4)
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_d_lapse', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_d_lapse for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_d_shift', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_d_shift for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_d_curv', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_d_curv for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_d_metric', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_d_metric for boundary condition' )
      end if
    end select
  
    select case ( sym_selector )
    case (1, 3)
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_res_lapse', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_res_lapse for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_res_shift', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_res_shift for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_res_curv', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_res_curv for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_res_metric', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_res_metric for boundary condition' )
      end if
    end select
  
    select case ( sym_selector )
    case (1, 2, 3)
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_red_all', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_red_all for boundary condition' )
      end if
    end select
  
    select case ( sym_selector )
    case (2)
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_q_lapse', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_q_lapse for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_q_shift', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_q_shift for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_q_curv', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_q_curv for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                             'NoExcision::cg_q_metric', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select cg_q_metric for boundary condition' )
      end if
    end select
  
    select case ( sym_selector )
    case (3)
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                                'ADMBase::lapse', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select lapse for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                                'ADMBase::shift', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select shift for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                                'ADMBase::curv', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select curv for boundary condition' )
      end if
  
      ierr = Boundary_SelectGroupForBC ( cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, &
                                                'ADMBase::metric', 'None' )
      if ( ierr /= 0 ) then
        call CCTK_WARN ( 0, 'Could not select metric for boundary condition' )
      end if
    end select

  end if

end subroutine NoExcision_CGApplySym


subroutine NoExcision_Set_Zero(CCTK_ARGUMENTS)

  use NoExcision_mod

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: cx, cy, cz, iradx, irady, iradz
  CCTK_REAL, dimension(:,:,:), allocatable :: dist2
  integer :: n
  CCTK_INT :: my_level, n_levels

  call NoExcision_levelinfo ( cctkGH, my_level, n_levels )

  if ( my_level == n_levels-1 ) then

    allocate ( dist2(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) )

    do n = 1, num_regions

      cx = centre_x(n)
      cy = centre_y(n)
      cz = centre_z(n)

      if (CCTK_EQUALS(region_shape(n), "sphere")) then

        iradx = 1 / radius(n)
        irady = 1 / radius(n)
        iradz = 1 / radius(n)

      else if (CCTK_EQUALS(region_shape(n), "ellipsoid")) then

        iradx = 1 / radius_x(n)
        irady = 1 / radius_y(n)
        iradz = 1 / radius_z(n)

      else

        call CCTK_WARN (0, "internal error")

      end if 

!$OMP PARALLEL WORKSHARE

      dist2 =   ((x - cx) * iradx)**2 + ((y - cy) * irady)**2 &
           &  + ((z - cz) * iradz)**2

      where ( dist2 <= one )

        gxx = zero
        gxy = zero
        gxz = zero
        gyy = zero
        gyz = zero
        gzz = zero
        kxx = zero
        kxy = zero
        kxz = zero
        kyy = zero
        kyz = zero
        kzz = zero
        alp = zero
        betax = zero
        betay = zero
        betaz = zero

      end where

!$OMP END PARALLEL WORKSHARE

    end do

    deallocate ( dist2 )

  end if
end subroutine NoExcision_Set_Zero
