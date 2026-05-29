! vim: syntax=fortran
#include "cctk.h"

module NullEvol_hyper_beta
contains

  subroutine NullEvol_beta (i, jns, bns, j_wt, b_wt, x_wt, mask)
    use NullGrid_Vars 
    implicit none

    CCTK_INT,                                   intent (in)    :: i
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: jns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (inout) :: bns 
    CCTK_INT,     dimension (lsh(1),lsh(2)),    intent (in)    :: mask
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (in)    :: j_wt
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (inout) :: b_wt, x_wt

    logical, save :: FirstTime = .true.

    CCTK_REAL,    dimension (:,:), allocatable, save::  K, dx_K
    CCTK_COMPLEX, dimension (:,:), allocatable, save::  J, dx_J

    CCTK_REAL xhere
    CCTK_INT i1, i2

    if (FirstTime) then
       FirstTime=.false.
       allocate(K(lsh(1),lsh(2)), dx_K(lsh(1),lsh(2)),&
            J(lsh(1),lsh(2)), dx_J(lsh(1),lsh(2)))
       K=0; dx_K=0; J=0; dx_J=0
    end if

    if(minval(mask).eq.0) then

       bns(:,:,i) = b_wt
       ! for points not to near the WT we improve this 
       J = 0.5 * (jns(:,:,i) + j_wt)
       K = sqrt(1. + J * conjg(J))

       do i2 = 1, lsh(2)
          do i1 = 1, lsh(1)
             if(abs(xb(i)-x_wt(i1,i2)).gt.1.e-5*dx) then

                ! on all points closer than dx to the boundary we use 
                ! a start-up algorithm with improved accuracy

                dx_J(i1,i2) = (jns(i1,i2,i) - j_wt(i1,i2))/(xb(i)-x_wt(i1,i2))
                dx_K(i1,i2) = dble( dx_J(i1,i2) * conjg(J(i1,i2)) ) / K(i1,i2)

                bns(i1,i2,i) = b_wt(i1,i2) + (xb(i)-x_wt(i1,i2)) * 0.5*(x_wt(i1,i2)+xb(i))&
                                                           * (1. - 0.5*(x_wt(i1,i2)+xb(i)))&
                          / 8. * (dx_J(i1,i2) * conjg(dx_J(i1,i2)) - dx_K(i1,i2) * dx_K(i1,i2))

             end if
          end do
       end do
    end if

    ! here comes the standard evolution algorithm

    xhere = xbh(i-1)

    J = 0.5 * (jns(:,:,i) + jns(:,:,i-1))
    K = sqrt(1. + J * conjg(J))

    dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx
    dx_K = dble( dx_J * conjg(J) ) / K

    bns(:,:,i) = bns(:,:,i) * (1 - mask) + mask * ( bns(:,:,i-1) &
         + dx * xhere * (1. - xhere) / 8. * (dx_J * conjg(dx_J) - dx_K * dx_K) )

  end subroutine NullEvol_beta

  subroutine NullEvol_cb (i, bns, cbns, mask)
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_REAL,    dimension (lsh(1), lsh(2), nx), intent (in)    :: bns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: cbns 
    CCTK_INT,     dimension (lsh(1), lsh(2)),     intent (in)    :: mask

    logical, save :: FirstTime = .true.

    CCTK_REAL,    dimension (:,:), allocatable, save ::  dx_beta
    CCTK_COMPLEX, dimension (:,:), allocatable, save ::  eth_dx_beta

    if (FirstTime) then
       FirstTime=.false.
       allocate(dx_beta(lsh(1),lsh(2)), eth_dx_beta(lsh(1),lsh(2)))
       dx_beta=0; eth_dx_beta=0
    end if

    if(minval(mask).eq.0) call NullEvol_Set_cb(i, bns, cbns)

    dx_beta = (bns(:,:,i) - bns(:,:,i-1)) / dx

    call NullInterp_d1(eth_dx_beta, dcmplx(dx_beta), 0_ik, 1_ik) 

    cbns(:,:,i) = cbns(:,:,i) * (1 - mask) &
         + mask * (cbns(:,:,i-1) + dx * eth_dx_beta)

  end subroutine NullEvol_cb

  subroutine NullEvol_ck (i, jns, ckns, mask)
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: ckns 
    CCTK_INT,     dimension (lsh(1),lsh(2)),      intent (in)    :: mask

    logical, save :: FirstTime = .true.
    CCTK_REAL,    dimension (:,:), allocatable, save :: K, dx_K
    CCTK_COMPLEX, dimension (:,:), allocatable, save :: J, dx_J, eth_dx_K

    if (FirstTime) then
       FirstTime=.false.
       allocate(K(lsh(1),lsh(2)), dx_K(lsh(1),lsh(2)), J(lsh(1),lsh(2)),&
            dx_J(lsh(1),lsh(2)), eth_dx_K(lsh(1),lsh(2)))
       K=0; dx_K=0; J=0; dx_J=0; eth_dx_K=0
    end if

    if(minval(mask).eq.0) call NullEvol_Set_ck(i, jns, ckns)

    J = 0.5 * (jns(:,:,i) + jns(:,:,i-1))
    K = sqrt(1. + J * conjg(J))

    dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx
    dx_K = dble( dx_J * conjg(J) ) / K

    call NullInterp_d1(eth_dx_K, dcmplx(dx_K), 0_ik, 1_ik) 

    ckns(:,:,i) = ckns(:,:,i) * (1 - mask) &
         + mask * (ckns(:,:,i-1) + dx * eth_dx_K)

  end subroutine NullEvol_ck

  subroutine NullEvol_nu (i, jns, nuns, mask)
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: nuns 
    CCTK_INT,     dimension (lsh(1), lsh(2)),     intent (in)    :: mask

    logical, save :: FirstTime = .true.
    CCTK_COMPLEX, dimension (:,:), allocatable, save ::  dx_J, ethb_dx_J

    if (FirstTime) then
       FirstTime=.false.
       allocate(dx_J(lsh(1),lsh(2)), ethb_dx_J(lsh(1),lsh(2)))
       dx_J=0; ethb_dx_J=0
    end if

    if(minval(mask).eq.0) call NullEvol_Set_nu(i, jns, nuns)

    dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx

    call NullInterp_d1(ethb_dx_J, dx_J, 2_ik, -1_ik) 

    nuns(:,:,i) = nuns(:,:,i) * (1 - mask)  + mask * (nuns(:,:,i-1) + dx * ethb_dx_J)

  end subroutine NullEvol_nu

  subroutine NullEvol_Set_cb (i, bns, cbns)
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_REAL,    dimension (lsh(1), lsh(2), nx), intent (in)    :: bns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: cbns 

    call NullInterp_d1(cbns(:,:,i), dcmplx(bns(:,:,i)), 0_ik, 1_ik) 

  end subroutine NullEvol_Set_cb


  subroutine NullEvol_Set_ck (i, jns, ckns)
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: ckns 

    CCTK_COMPLEX, dimension (lsh(1), lsh(2)) :: K

    K = sqrt(1. + jns(:,:,i) * conjg(jns(:,:,i)))
    call NullInterp_d1(ckns(:,:,i), K, 0_ik, 1_ik) 

  end subroutine NullEvol_Set_ck

  subroutine NullEvol_Set_nu (i, jns, nuns)
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: nuns 

    call NullInterp_d1(nuns(:,:,i), jns(:,:,i), 2_ik, -1_ik) 

  end subroutine NullEvol_Set_nu

end module NullEvol_hyper_beta
