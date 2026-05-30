! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

module NullInterp_Deriv
use cctk
implicit none
contains

   subroutine NullInterp_wt_da_ethr(cctkGH, tmp_cgfn, tmp_cgfs, fn, fs, fn_q, fs_q, fn_p, fs_p)

   use NullGrid_Vars
   use NullInterp_Interp
   use NullInterp_Eth
    implicit none

    CCTK_POINTER, intent(in) :: cctkGH
    CCTK_REAL, dimension (lsh(1),lsh(2)), intent (in)  :: fn, fs
    CCTK_REAL, dimension (lsh(1),lsh(2)), intent (out) :: fn_q, fs_q, fn_p, fs_p
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)) :: eth_fn, eth_fs
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)) :: ethb_fn, ethb_fs
    CCTK_COMPLEX, dimension(lsh(1),lsh(2)), intent(inout) :: tmp_cgfn, tmp_cgfs
    CCTK_COMPLEX :: i
    CCTK_INT :: zero = 0
    CCTK_INT :: one = 1
    CCTK_INT :: m1 = -1

    i = dcmplx(0.,1.)
    call NullInterp_d1 (eth_fn, dcmplx(fn), zero, one)
    call NullInterp_d1 (eth_fs, dcmplx(fs), zero, one)
    call NullInterp_d1 (ethb_fn, dcmplx(fn), zero, m1)
    call NullInterp_d1 (ethb_fs, dcmplx(fs), zero, m1)
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_fn, eth_fs, one)
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, ethb_fn, ethb_fs, m1)
  
    fn_q(1:lsh(1),1:lsh(2)) =&
        dble(eth_fn(1:lsh(1),1:lsh(2)) +&
         ethb_fn(1:lsh(1),1:lsh(2))) / pp
    fs_q(1:lsh(1),1:lsh(2)) =&
        dble(eth_fs(1:lsh(1),1:lsh(2)) +&
        ethb_fs(1:lsh(1),1:lsh(2))) / pp

    fn_p(1:lsh(1),1:lsh(2)) =&
         dble(i * (ethb_fn(1:lsh(1),1:lsh(2)) -&
         eth_fn(1:lsh(1),1:lsh(2)))) / pp
    fs_p(1:lsh(1),1:lsh(2)) =&
         dble(i * (ethb_fs(1:lsh(1),1:lsh(2)) -&
          eth_fs(1:lsh(1),1:lsh(2)))) / pp

end subroutine NullInterp_wt_da_ethr

subroutine NullInterp_wt_d2a_ethr(cctkGH, tmp_cgfn, tmp_cgfs,&
                                 fn, fs, fn_q, fs_q, fn_p, fs_p,&
                                 fn_qq, fs_qq, fn_qp, fs_qp, fn_pp, fs_pp)
use NullGrid_Vars
use NullInterp_Interp
use NullInterp_Eth

    implicit none

    CCTK_POINTER, intent(in) :: cctkGH
    CCTK_REAL, dimension (lsh(1),lsh(2)), intent (in)  :: fn, fs
    CCTK_REAL, dimension (lsh(1),lsh(2)), intent (out) :: fn_q, fs_q, fn_p, fs_p,&
                                          fn_qq, fs_qq, fn_qp, fs_qp, fn_pp, fs_pp
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)) :: eth_fn, eth_fs,&
                             ethb_fn, ethb_fs, eth_eth_fn, eth_eth_fs, &
                             eth_ethb_fn, eth_ethb_fs, ethb_ethb_fn, ethb_ethb_fs
    CCTK_COMPLEX, dimension(lsh(1),lsh(2)), intent(inout) :: tmp_cgfn, tmp_cgfs    
    CCTK_INT :: zero = 0
    CCTK_INT :: one = 1
    CCTK_INT :: two = 2
    CCTK_INT :: m1 = -1
    CCTK_INT :: m2 = -2

    CCTK_COMPLEX :: ii = (0., 1.)

    call NullInterp_d1 (eth_fn, dcmplx(fn), zero, one)
    call NullInterp_d1 (eth_fs, dcmplx(fs), zero, one)
    call NullInterp_d1 (ethb_fn, dcmplx(fn), zero, m1)
    call NullInterp_d1 (ethb_fs, dcmplx(fs), zero, m1)
    
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_fn, eth_fs, one)
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_fn, eth_fs, m1)

    call NullInterp_d2 (eth_eth_fn,   dcmplx(fn), zero,  one,  one)
    call NullInterp_d2 (eth_eth_fs,   dcmplx(fs), zero,  one,  one)
    call NullInterp_d2 (eth_ethb_fn,  dcmplx(fn), zero,  one, m1)
    call NullInterp_d2 (eth_ethb_fs,  dcmplx(fs), zero,  one, m1)
    call NullInterp_d2 (ethb_ethb_fn, dcmplx(fn), zero, m1, m1)
    call NullInterp_d2 (ethb_ethb_fs, dcmplx(fs), zero, m1, m1)
 
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_eth_fn, eth_eth_fs, two)
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_ethb_fn, eth_ethb_fs, zero)
    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, ethb_ethb_fn, ethb_ethb_fs, m2)

    fn_q(1:lsh(1),1:lsh(2)) =&
       dble(eth_fn(1:lsh(1),1:lsh(2)) +&
       ethb_fn(1:lsh(1),1:lsh(2))) / pp
    fs_q(1:lsh(1),1:lsh(2)) =&
       dble(eth_fs(1:lsh(1),1:lsh(2)) +&
       ethb_fs(1:lsh(1),1:lsh(2))) / pp

    fn_p(1:lsh(1),1:lsh(2)) =&
       dble(ii * (ethb_fn(1:lsh(1),1:lsh(2)) -&
       eth_fn(1:lsh(1),1:lsh(2)))) / pp
    fs_p(1:lsh(1),1:lsh(2)) =&
       dble(ii * (ethb_fs(1:lsh(1),1:lsh(2)) -&
       eth_fs(1:lsh(1),1:lsh(2)))) / pp

    fn_qq(1:lsh(1),1:lsh(2)) =&
       dble(-(-ethb_ethb_fn(1:lsh(1),1:lsh(2))+&
       2*qs*ethb_fn(1:lsh(1),1:lsh(2))-&
       2*ii*ps*ethb_fn(1:lsh(1),1:lsh(2))-&
       2*eth_ethb_fn(1:lsh(1),1:lsh(2))-&
       eth_eth_fn(1:lsh(1),1:lsh(2))+&
       2*qs*eth_fn(1:lsh(1),1:lsh(2))+&
       2*ii*ps*eth_fn(1:lsh(1),1:lsh(2)))/pp**2)

    fs_qq(1:lsh(1),1:lsh(2)) =&
       dble(-(-ethb_ethb_fs(1:lsh(1),1:lsh(2))+&
       2*qs*ethb_fs(1:lsh(1),1:lsh(2))-&
       2*ii*ps*ethb_fs(1:lsh(1),1:lsh(2))-&
       2*eth_ethb_fs(1:lsh(1),1:lsh(2))-&
       eth_eth_fs(1:lsh(1),1:lsh(2))+&
       2*qs*eth_fs(1:lsh(1),1:lsh(2))+&
       2*ii*ps*eth_fs(1:lsh(1),1:lsh(2)))/pp**2)

    fn_qp(1:lsh(1),1:lsh(2)) =&
       dble((ii*ethb_ethb_fn(1:lsh(1),1:lsh(2))-&
       2*ii*qs*ethb_fn(1:lsh(1),1:lsh(2))-&
       2*ps*ethb_fn(1:lsh(1),1:lsh(2))-&
       ii*eth_eth_fn(1:lsh(1),1:lsh(2))+&
       2*ii*qs*eth_fn(1:lsh(1),1:lsh(2))-&
       2*ps*eth_fn(1:lsh(1),1:lsh(2)))/pp**2)

    fs_qp(1:lsh(1),1:lsh(2)) =&
       dble((ii*ethb_ethb_fs(1:lsh(1),1:lsh(2))-&
       2*ii*qs*ethb_fs(1:lsh(1),1:lsh(2))-&
       2*ps*ethb_fs(1:lsh(1),1:lsh(2))-&
       ii*eth_eth_fs(1:lsh(1),1:lsh(2))+&
       2*ii*qs*eth_fs(1:lsh(1),1:lsh(2))-&
       2*ps*eth_fs(1:lsh(1),1:lsh(2)))/pp**2)

    fn_pp(1:lsh(1),1:lsh(2)) =&
       dble((2*eth_ethb_fn(1:lsh(1),1:lsh(2))-&
       ethb_ethb_fn(1:lsh(1),1:lsh(2))+&
       2*qs*ethb_fn(1:lsh(1),1:lsh(2))-&
       2*ii*ps*ethb_fn(1:lsh(1),1:lsh(2))-&
       eth_eth_fn(1:lsh(1),1:lsh(2))+&
       2*qs*eth_fn(1:lsh(1),1:lsh(2))+&
       2*ii*ps*eth_fn(1:lsh(1),1:lsh(2)))/pp**2)

    fs_pp(1:lsh(1),1:lsh(2)) =&
       dble((2*eth_ethb_fs(1:lsh(1),1:lsh(2))-&
       ethb_ethb_fs(1:lsh(1),1:lsh(2))+&
       2*qs*ethb_fs(1:lsh(1),1:lsh(2))-&
       2*ii*ps*ethb_fs(1:lsh(1),1:lsh(2))-&
       eth_eth_fs(1:lsh(1),1:lsh(2))+&
       2*qs*eth_fs(1:lsh(1),1:lsh(2))+&
       2*ii*ps*eth_fs(1:lsh(1),1:lsh(2)))/pp**2)
end subroutine NullInterp_wt_d2a_ethr

subroutine NullInterp_wt_da(nq, np, dd, f, fq, fp)
!use NullInterp_Interp
!use NullInterp_Eth

   implicit none

   CCTK_INT,                    intent(in)  :: nq, np
   CCTK_REAL,                   intent(in)  :: dd
   CCTK_REAL, dimension(nq,np), intent(in)  :: f
   CCTK_REAL, dimension(nq,np), intent(out) :: fq, fp

   CCTK_INT :: i, j
   CCTK_REAL :: h1

! assume dq = dp :(
   h1 = dble(1. / (2. * dd))

   do j = 1, np

      do i = 2, nq-1
         fq(i,j) = dble((f(i+1,j) - f(i-1,j)) * h1)
      end do
      fq(1,j)  = dble((-3. * f(1,j) + 4. * f(2,j) - f(3,j)) * h1)
      fq(nq,j) = dble(( 3. * f(nq,j) - 4. * f(nq-1,j) + f(nq-2,j)) * h1)

   end do

   do i = 1, nq

      do j = 2, np-1
         fp(i,j) = dble((f(i,j+1) - f(i,j-1)) * h1)
      end do
      fp(i,1)  = dble((-3. * f(i,1) + 4. * f(i,2) - f(i,3)) * h1)
      fp(i,np) = dble(( 3. * f(i,np) - 4. * f(i,np-1) + f(i,np-2)) * h1)

   end do

end subroutine NullInterp_wt_da

end module NullInterp_Deriv
