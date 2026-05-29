! vim: syntax=fortran

#include "cctk.h"

module NullConstr_Util

contains

!define DEBUG_I 4

subroutine x2r_derivs(dr_dx,d2r_dx2,dr_dxh,d2r_dx2h,null_nx,xb,xbh)
  use NullGrid_Vars, only: rwt
  implicit none
  CCTK_INT, intent(in) :: null_nx
  CCTK_REAL, dimension(null_nx), intent(in) :: xb, xbh
  CCTK_REAL, dimension(null_nx), intent(out) :: dr_dx,d2r_dx2,dr_dxh,d2r_dx2h
  
  ! For r=rwt*x/(1-x)

  dr_dx=rwt/(1-xb)**2
  d2r_dx2=2*rwt/(1-xb)**3

  dr_dxh=rwt/(1-xbh)**2
  d2r_dx2h=2*rwt/(1-xbh)**3
 
  ! FIXME: this is only first order accurate
 
  dr_dx(null_nx)=rwt/(1-xb(null_nx-1))**2
  d2r_dx2(null_nx)=2*rwt/(1-xb(null_nx-1))**3

  dr_dxh(null_nx)=rwt/(1-xbh(null_nx-1))**2
  d2r_dx2h(null_nx)=2*rwt/(1-xbh(null_nx-1))**3
  
end subroutine x2r_derivs

subroutine real_derivs (aa_n,aa_p,dt,dx,null_nx,null_lsh,i,dr_dx,d2r_dx2, &
  aa_00,aa_01,aa_02,aa_03,aa_04,aa_11,aa_12,aa_13,aa_14, &
  aa_22,aa_23,aa_24,aa_33,aa_34)
  use NullInterp
  implicit none
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_INT, intent(in)  :: null_nx,null_lsh(2),i
  CCTK_REAL, intent(in) :: dt,dx,dr_dx,d2r_dx2
  CCTK_REAL, dimension(null_lsh(1),null_lsh(2),null_nx), intent(in) :: aa_n,aa_p
  CCTK_REAL, dimension(null_lsh(1),null_lsh(2)), intent(out) :: aa_00,aa_01,aa_04,aa_11,aa_14
  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2)), intent(out)  :: aa_02,aa_03,aa_12,aa_13,aa_22, &  ! changed (added): intent(out)
  aa_23,aa_24,aa_33,aa_34
  CCTK_REAL, dimension(null_lsh(1),null_lsh(2)) :: aa_0x,aa_xx,aa_x4
 
  if (i.le.1.or.i.gt.null_nx) call CCTK_WARN(0, "invalid index")

  aa_00=(aa_n(:,:,i)+aa_p(:,:,i))/2
  aa_04=(aa_n(:,:,i)-aa_p(:,:,i))/dt

  if (i.eq.null_nx) then ! FIXME: for now use only first order accuracy at scri
     aa_0x=(aa_n(:,:,i)-aa_n(:,:,i-1)+aa_p(:,:,i)-aa_p(:,:,i-1))/(2*dx)
     aa_01=aa_0x/dr_dx
     aa_x4=(aa_n(:,:,i)-aa_n(:,:,i-1)-aa_p(:,:,i)+aa_p(:,:,i-1))/(dt*dx)
     aa_14=aa_x4/dr_dx
     aa_xx=(aa_n(:,:,(i-1)+1)-2*aa_n(:,:,(i-1))+aa_n(:,:,(i-1)-1) &
     +aa_p(:,:,(i-1)+1)-2*aa_p(:,:,(i-1))+aa_p(:,:,(i-1)-1))/(2*dx**2)
     aa_11=aa_xx/dr_dx**2-aa_0x*d2r_dx2/dr_dx**3
  else
     aa_0x=(aa_n(:,:,i+1)-aa_n(:,:,i-1)+aa_p(:,:,i+1)-aa_p(:,:,i-1))/(4*dx)
     aa_01=aa_0x/dr_dx
     aa_x4=(aa_n(:,:,i+1)-aa_n(:,:,i-1)-aa_p(:,:,i+1)+aa_p(:,:,i-1))/(2*dt*dx)
     aa_14=aa_x4/dr_dx
     aa_xx=(aa_n(:,:,i+1)-2*aa_n(:,:,i)+aa_n(:,:,i-1) &
     +aa_p(:,:,i+1)-2*aa_p(:,:,i)+aa_p(:,:,i-1))/(2*dx**2)
     aa_11=aa_xx/dr_dx**2-aa_0x*d2r_dx2/dr_dx**3
  end if

  call NullInterp_d1(aa_02, dcmplx(aa_00,0), 0_ik, 1_ik) 
  call NullInterp_d1(aa_03, dcmplx(aa_00,0), 0_ik, -1_ik) 
  call NullInterp_d1(aa_12, dcmplx(aa_01,0), 0_ik, 1_ik) 
  call NullInterp_d1(aa_13, dcmplx(aa_01,0), 0_ik, -1_ik) 
  call NullInterp_d1(aa_24, dcmplx(aa_04,0), 0_ik, 1_ik) 
  call NullInterp_d1(aa_34, dcmplx(aa_04,0), 0_ik, -1_ik) 
  call NullInterp_d1(aa_22, aa_02, 1_ik, 1_ik) 
  call NullInterp_d1(aa_33, aa_03, -1_ik, -1_ik) 
  call NullInterp_d1(aa_23, aa_03, -1_ik, 1_ik) 

#ifdef DEBUG_I
  if (i.eq.DEBUG_I) then 
    write (*,*) 'real input:', maxval(abs(aa_n(:,:,i))), maxval(abs(aa_p(:,:,i)))
    write (*,*) 'real output1:', maxval(abs(aa_00)),maxval(abs(aa_01)),maxval(abs(aa_04)),maxval(abs(aa_11)),maxval(abs(aa_14))
    write (*,*) 'real output2:', maxval(abs(aa_02)),maxval(abs(aa_03)),maxval(abs(aa_12)),maxval(abs(aa_13)),&
      maxval(abs(aa_22)),maxval(abs(aa_23)),maxval(abs(aa_24)),maxval(abs(aa_33)),maxval(abs(aa_34))
  end if
#endif

end subroutine real_derivs

subroutine cmplx_derivs &
  (aa_n,aa_p,dt,dx,null_nx,null_lsh,i,spin_wt,dr_dx,d2r_dx2,half_grid, &
  aa_00,aa_01,aa_02,aa_03,aa_04,aa_11,aa_12,aa_13,aa_14, &
  aa_22,aa_23,aa_24,aa_33,aa_34, &
  aab_00,aab_01,aab_02,aab_03,aab_04,aab_11,aab_12,aab_13,aab_14, &
  aab_22,aab_23,aab_24,aab_33,aab_34)
  use NullInterp
  implicit none
  
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_INT, intent(in)  :: null_nx,null_lsh(2),i,spin_wt
  LOGICAL, intent(in) :: half_grid
  CCTK_REAL, intent(in) :: dt,dx,dr_dx,d2r_dx2

  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2),null_nx), intent(in) :: aa_n,aa_p
  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2)), intent(out) :: aa_00,aa_01,aa_04,aa_11,aa_14
  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2)), intent(out) :: aa_02,aa_03,aa_12,aa_13,aa_22, &
  aa_23,aa_24,aa_33,aa_34
  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2)), intent(out) :: aab_00,aab_01,aab_04,aab_11,aab_14
  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2)), intent(out) :: aab_02,aab_03,aab_12,aab_13,aab_22, &
  aab_23,aab_24,aab_33,aab_34
  CCTK_COMPLEX, dimension(null_lsh(1),null_lsh(2)):: aa_0x,aa_x4,aa_xx

  if (i.le.1.or.i.gt.null_nx) call CCTK_WARN(0, "invalid index")

  if(half_grid) then
     aa_00=(aa_n(:,:,i)+aa_p(:,:,i)+aa_n(:,:,i-1)+aa_p(:,:,i-1))/4
     aa_04=(aa_n(:,:,i)-aa_p(:,:,i)+aa_n(:,:,i-1)-aa_p(:,:,i-1))/(2*dt)
     aa_0x=(aa_n(:,:,i)-aa_n(:,:,i-1)+aa_p(:,:,i)-aa_p(:,:,i-1))/(2*dx)
     aa_x4=(aa_n(:,:,i)-aa_n(:,:,i-1)- &
     aa_p(:,:,i)+aa_p(:,:,i-1))/(dt*dx)
     if (i.eq.null_nx) then ! FIXME: for now only 1st order accurate
        aa_xx=(aa_n(:,:,i-1+1)-aa_n(:,:,i-1)-aa_n(:,:,i-1-1)+aa_n(:,:,i-1-2) &
        +aa_p(:,:,i-1+1)-aa_p(:,:,i-1)-aa_p(:,:,i-1-1)+aa_p(:,:,i-1-2) )/(4*dx**2)
     else
        aa_xx=(aa_n(:,:,i+1)-aa_n(:,:,i)-aa_n(:,:,i-1)+aa_n(:,:,i-2) &
        +aa_p(:,:,i+1)-aa_p(:,:,i)-aa_p(:,:,i-1)+aa_p(:,:,i-2) )/(4*dx**2)
     end if
  else
     aa_00=(aa_n(:,:,i)+aa_p(:,:,i))/2
     aa_04=(aa_n(:,:,i)-aa_p(:,:,i))/dt
     if (i.eq.null_nx) then  ! FIXME: for now only 1st order accurate
        aa_0x=(aa_n(:,:,i)-aa_n(:,:,i-1)+aa_p(:,:,i)-aa_p(:,:,i-1))/(2*dx)
        aa_x4=(aa_n(:,:,i)-aa_n(:,:,i-1)- &
        aa_p(:,:,i)+aa_p(:,:,i-1))/(dt*dx)
        aa_xx=(aa_n(:,:,i-1+1)-2*aa_n(:,:,i-1)+aa_n(:,:,i-1-1) &
        +aa_p(:,:,i-1+1)-2*aa_p(:,:,i-1)+aa_p(:,:,i-1-1))/(2*dx**2)
     else
        aa_0x=(aa_n(:,:,i+1)-aa_n(:,:,i-1)+aa_p(:,:,i+1)-aa_p(:,:,i-1))/(4*dx)
        aa_x4=(aa_n(:,:,i+1)-aa_n(:,:,i-1)- &
        aa_p(:,:,i+1)+aa_p(:,:,i-1))/(2*dt*dx)
        aa_xx=(aa_n(:,:,i+1)-2*aa_n(:,:,i)+aa_n(:,:,i-1) &
        +aa_p(:,:,i+1)-2*aa_p(:,:,i)+aa_p(:,:,i-1))/(2*dx**2)
     end if
  end if
  aa_01=aa_0x/dr_dx
  aa_14=aa_x4/dr_dx
  aa_11=aa_xx/dr_dx**2-aa_0x*d2r_dx2/dr_dx**3
  call NullInterp_d1(aa_02, aa_00, spin_wt, 1_ik) 
  call NullInterp_d1(aa_03, aa_00, spin_wt, -1_ik) 
  call NullInterp_d1(aa_12, aa_01, spin_wt, 1_ik) 
  call NullInterp_d1(aa_13, aa_01, spin_wt, -1_ik) 
  call NullInterp_d1(aa_24, aa_04, spin_wt, 1_ik) 
  call NullInterp_d1(aa_34, aa_04, spin_wt, -1_ik) 
  call NullInterp_d1(aa_22, aa_02, spin_wt+1, 1_ik) 
  call NullInterp_d1(aa_33, aa_03, spin_wt-1, -1_ik) 
  call NullInterp_d1(aa_23, aa_03, spin_wt-1, 1_ik)
  
  aab_00=conjg(aa_00)
  aab_01=conjg(aa_01)
  aab_04=conjg(aa_04)
  aab_14=conjg(aa_14)
  aab_11=conjg(aa_11)
  aab_02=conjg(aa_03)
  aab_03=conjg(aa_02)
  aab_12=conjg(aa_13)
  aab_13=conjg(aa_12)
  aab_24=conjg(aa_34)
  aab_34=conjg(aa_24)
  aab_22=conjg(aa_33)
  aab_33=conjg(aa_22)
  call NullInterp_d1(aab_23, aab_03, -spin_wt-1, 1_ik)

#ifdef DEBUG_I
  if (i.eq.DEBUG_I) then 
     write (*,*) 'cmplex input:', maxval(abs(aa_n(:,:,i))), maxval(abs(aa_p(:,:,i)))
     write (*,*) 'cmplex output3:', maxval(abs(aa_00)),maxval(abs(aa_01)),maxval(abs(aa_04)),maxval(abs(aa_11)),maxval(abs(aa_14))
     write (*,*) 'cmplex output4:', maxval(abs(aa_02)),maxval(abs(aa_03)),maxval(abs(aa_12)),maxval(abs(aa_13)),&
       maxval(abs(aa_22)),maxval(abs(aa_23)),maxval(abs(aa_24)),maxval(abs(aa_33)),maxval(abs(aa_34))
  end if
#endif

end subroutine cmplx_derivs


!subroutine NullInterp_d1(F_out, F_in, spin, e1)
!   implicit none
!   CCTK_COMPLEX, dimension (:,:), intent (in) :: F_in
!   CCTK_COMPLEX, dimension (:,:), intent (out) :: F_out
!   CCTK_INT, intent (in)  :: spin, e1

!  F_out=dcmplx(0.0,0.0)

!end subroutine NullInterp_d1

end module NullConstr_Util
