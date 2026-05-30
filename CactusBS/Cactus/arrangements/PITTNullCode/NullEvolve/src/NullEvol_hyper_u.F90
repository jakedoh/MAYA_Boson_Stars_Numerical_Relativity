! vim: syntax=fortran
#include "cctk.h"

module NullEvol_hyper_u
contains
  subroutine NullEvol_u (i, jns, bns, qns, uns, u_wt, u_x_wt, x_wt, mask)
    use NullEvol_Mask
    use NullInterp
    use NullGrid_Vars 
    implicit none

    CCTK_INT,                                   intent (in)    :: i
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: jns, qns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (in)    :: bns 
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (inout) :: uns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (in)    :: u_wt, u_x_wt
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: x_wt
    CCTK_INT,     dimension (lsh(1),lsh(2)),    intent (in)    :: mask

    CCTK_COMPLEX ::  x2Uc_x, uE, x2UE_x
    CCTK_REAL    ::  xc, xE
    CCTK_INT     :: i1, i2

    ! these store x^2 * U_{,x}
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)) :: x2U_x

    x2U_x = exp(2. * bns(:,:,i)) / rwt * (&
                                - jns(:,:,i) * conjg(qns(:,:,i))& ! J * conj(Q)
        + sqrt(1. + jns(:,:,i) * conjg(jns(:,:,i))) * qns(:,:,i)) ! K * Q

    if(minval(mask).eq.0) then

       do i2 = 1, lsh(2)
          do i1 = 1, lsh(1)
             ! We only change points between the WT and x(B+1/2).
             ! The points inside the WT are filled via a linear
             ! expansion in x from within the extraction algorithm.

             if(mask(i1,i2).eq.0) then

                xE = x_wt(i1,i2)
                uE = u_wt(i1,i2)
                x2UE_x = x_wt(i1,i2)**2 * u_x_wt(i1,i2)

                if(abs(xbh(i-1)-xE).lt.1.e-5*dx) then
                   uns(i1,i2,i-1) = u_wt(i1,i2)
                else if(xE < xbh(i-1)) then
                   ! we update the points i-1/2 based on dx_U(xE), xE and dx_U(i)
                   ! then we update i+1/2 based on i-1/2 and dx_U(i)
   
                   ! x(i-1/2) is outside the WT -- we update it
                   xc = (xE + xbh(i-1))/2
                   x2Uc_x = x2UE_x * (xb(i) -xc)/(xb(i)-xE) &
                          + x2U_x(i1,i2)  * (xE-xc)/(xE-xb(i))
                   uns(i1,i2,i-1) = uE + x2Uc_x * ( xbh(i-1)-xE ) / (xbh(i-1)*xE)
                end if

                if(abs(xbh(i)-xE).lt.1.e-5*dx) then
                   uns(i1,i2,i) = u_wt(i1,i2)
                else
                   ! we update the points i+1/2 based on dx_U(xE), xE and dx_U(i)
                   ! x(i-1/2) is outside the WT -- we update it
                   xc = (xE + xbh(i))/2
                   x2Uc_x = x2UE_x * (xb(i)-xc)/(xb(i)-xE)&
                          + x2U_x(i1,i2)  * (xE-xc)/(xE-xb(i))
                   uns(i1,i2,i) = uE + x2Uc_x * ( xbh(i)-xE ) / (xbh(i)*xE)
                end if

             end if ! mask

          end do
       end do
    end if

    ! update U(x(i+1/2))
    uns(:,:,i) = uns(:,:,i) * (1 - mask) + mask * (uns(:,:,i-1) + dx*x2U_x/xbh(i-1)/xbh(i))

  end subroutine NullEvol_u
end module NullEvol_hyper_u   
