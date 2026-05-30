#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"



subroutine LWT_calc_rhs (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
   CCTK_INT, parameter :: izero = 0
   CCTK_INT :: i, j, k
   CCTK_REAL :: rpt

   integer, parameter :: ik = kind (izero)
   
   integer :: nx, ny, nz
   CCTK_REAL, dimension(:,:,:), allocatable &
        :: hxx, hxy, hxz, hyy, hyz, hzz, &
           vux, vuy, vuz, &
           dxvux, dyvuy, dzvuz, &
           davux, dbvux, dcvux, &
           davuy, dbvuy, dcvuy, &
           davuz, dbvuz, dcvuz, &
           dxrho, dyrho, dzrho, &
           dxu, dyu, dzu, &
           dxxu, dxyu, dxzu, dyyu, dyzu, dzzu, &
           darho, dbrho, dcrho, &
           wuxx, wuxy, wuxz, wuyy, wuyz, wuzz, &
           dxwuxx, dxwuxy, dxwuxz, dxwuyy, dxwuyz, dxwuzz, &
           dywuxx, dywuxy, dywuxz, dywuyy, dywuyz, dywuzz, &
           dzwuxx, dzwuxy, dzwuxz, dzwuyy, dzwuyz, dzwuzz
 
   nx = cctk_lsh(1); ny = cctk_lsh(2); nz = cctk_lsh(3)
 
   allocate ( hxx(nx,ny,nz), hxy(nx,ny,nz), hxz(nx,ny,nz), &
              hyy(nx,ny,nz), hyz(nx,ny,nz), hzz(nx,ny,nz), &
              vux(nx,ny,nz), vuy(nx,ny,nz), vuz(nx,ny,nz), &
              wuxx(nx,ny,nz), wuxy(nx,ny,nz), wuxz(nx,ny,nz), &
              wuyy(nx,ny,nz), wuyz(nx,ny,nz), wuzz(nx,ny,nz), &
              dxwuxx(nx,ny,nz), dxwuxy(nx,ny,nz), dxwuxz(nx,ny,nz), &
              dxwuyy(nx,ny,nz), dxwuyz(nx,ny,nz), dxwuzz(nx,ny,nz), &
              dywuxx(nx,ny,nz), dywuxy(nx,ny,nz), dywuxz(nx,ny,nz), &
              dywuyy(nx,ny,nz), dywuyz(nx,ny,nz), dywuzz(nx,ny,nz), &
              dzwuxx(nx,ny,nz), dzwuxy(nx,ny,nz), dzwuxz(nx,ny,nz), &
              dzwuyy(nx,ny,nz), dzwuyz(nx,ny,nz), dzwuzz(nx,ny,nz) &
            )
   
   ! H^ij = g^ij alpha^2 - beta^i beta^j
   ! TODO: maybe change H^ij -> H^ij / alpha**2 for consistency
   

! NOTE: omp parallel workshare is commented out in this routine as it
! leads to wrong results on ifort

!!!$OMP PARALLEL WORKSHARE
   hxx = guxx * alpha**2 - betax * betax
   hxy = guxy * alpha**2 - betax * betay
   hxz = guxz * alpha**2 - betax * betaz
   hyy = guyy * alpha**2 - betay * betay
   hyz = guyz * alpha**2 - betay * betaz
   hzz = guzz * alpha**2 - betaz * betaz
   
   ! vu^i = epsilon / alpha (beta^i rho + H^ij v_j)
   vux = epsilon / alpha * (betax*rho)
   vuy = epsilon / alpha * (betay*rho)
   vuz = epsilon / alpha * (betaz*rho)
 
   wuxx = epsilon / alpha * hxx
   wuxy = epsilon / alpha * hxy
   wuxz = epsilon / alpha * hxz
   wuyy = epsilon / alpha * hyy
   wuyz = epsilon / alpha * hyz
   wuzz = epsilon / alpha * hzz
!!!$OMP END PARALLEL WORKSHARE
   
   deallocate ( hxx, hxy, hxz, hyy, hyz, hzz )
 
   allocate ( dxvux(nx,ny,nz), dyvuy(nx,ny,nz), dzvuz(nx,ny,nz) )
   allocate ( davux(nx,ny,nz), davuy(nx,ny,nz), davuz(nx,ny,nz) )
   allocate ( dbvux(nx,ny,nz), dbvuy(nx,ny,nz), dbvuz(nx,ny,nz) )
   allocate ( dcvux(nx,ny,nz), dcvuy(nx,ny,nz), dcvuz(nx,ny,nz) )
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Derivatives are done here. Substitute the correct Jacobians
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call globalDiff_gv (cctkGH, 0_ik, vux, dxvux, J11, J21, J31, &
                                  J12, J22, J32, J13, J23, J33,-1_ik)
   call globalDiff_gv (cctkGH, 1_ik, vuy, dyvuy, J11, J21, J31, &
                                  J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 2_ik, vuz, dzvuz, J11, J21, J31, &
                                  J12, J22, J32, J13, J23, J33, -1_ik)
 
 
   call globalDiff_gv (cctkGH, 0_ik, wuxx, dxwuxx, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 0_ik, wuxy, dxwuxy, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 0_ik, wuxz, dxwuxz, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)

   call globalDiff_gv (cctkGH, 1_ik, wuxy, dywuxy, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 1_ik, wuyy, dywuyy, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 1_ik, wuyz, dywuyz, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
 
   call globalDiff_gv (cctkGH, 2_ik, wuxz, dzwuxz, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 2_ik, wuyz, dzwuyz, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 2_ik, wuzz, dzwuzz, J11, J21, J31, &
                                    J12, J22, J32, J13, J23, J33, -1_ik)
 
   deallocate ( vux, vuy, vuz )
        
   allocate ( dxrho(nx,ny,nz), dyrho(nx,ny,nz), dzrho(nx,ny,nz) )
   allocate ( dxu(nx,ny,nz), dyu(nx,ny,nz), dzu(nx,ny,nz) )
   allocate ( dxxu(nx,ny,nz), dxyu(nx,ny,nz), dxzu(nx,ny,nz), &
              dyyu(nx,ny,nz), dyzu(nx,ny,nz), dzzu(nx,ny,nz) )
   allocate ( darho(nx,ny,nz), dbrho(nx,ny,nz), dcrho(nx,ny,nz) )

!!!!!!!!!!!!!!!!!!! More Derivatives. Substitute Jacobians again !!
   call globalDiff_gv (cctkGH, 0_ik, rho, dxrho, J11, J21, J31, &
                                  J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 1_ik, rho, dyrho, J11, J21, J31, &
                                  J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 2_ik, rho, dzrho, J11, J21, J31, &
                                  J12, J22, J32, J13, J23, J33, -1_ik)
 
!!!$OMP PARALLEL WORKSHARE
   dx_rho = dxrho
   dy_rho = dyrho
   dz_rho = dzrho
!!!$OMP END PARALLEL WORKSHARE
 
   call globalDiff_gv (cctkGH, 0_ik, u, dxu, J11, J21, J31, &
                              J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 1_ik, u, dyu, J11, J21, J31, &
                              J12, J22, J32, J13, J23, J33, -1_ik)
   call globalDiff_gv (cctkGH, 2_ik, u, dzu, J11, J21, J31, &
                              J12, J22, J32, J13, J23, J33, -1_ik)

!!!$OMP PARALLEL WORKSHARE
   dx_u = dxu
   dy_u = dyu
   dz_u = dzu
!!!$OMP END PARALLEL WORKSHARE

   if (compute_second_derivative_from_first_derivative/=0) then

     write(*,*) "WARNING: Computing second derivates from first derivatives"

     call globalDiff_gv (cctkGH, 0_ik, dxu, dxxu, J11, J21, J31, &
                                 J12, J22, J32, J13, J23, J33, -1_ik)
     call globalDiff_gv (cctkGH, 1_ik, dxu, dxyu, J11, J21, J31, &
                                   J12, J22, J32, J13, J23, J33, -1_ik)
     call globalDiff_gv (cctkGH, 2_ik, dxu, dxzu, J11, J21, J31, &
                                   J12, J22, J32, J13, J23, J33, -1_ik)
     call globalDiff_gv (cctkGH, 1_ik, dyu, dyyu, J11, J21, J31, &
                                   J12, J22, J32, J13, J23, J33, -1_ik)
     call globalDiff_gv (cctkGH, 2_ik, dyu, dyzu, J11, J21, J31, &
                                   J12, J22, J32, J13, J23, J33, -1_ik)
     call globalDiff_gv (cctkGH, 2_ik, dzu, dzzu, J11, J21, J31, &
                                   J12, J22, J32, J13, J23, J33, -1_ik)

   else

!!!!!!!!!!!!!!!!! Subs Jacobians   !!!
     call globalDiff2_gv (cctkGH, 0_ik, 0_ik, u, dxxu, J11, J21, J31, &
                                        J12, J22, J32, J13, J23, J33, &
                                        dJ111, dJ211, dJ311, &
                                        dJ112, dJ212, dJ312, &
                                        dJ113, dJ213, dJ313, &
                                        dJ122, dJ222, dJ322, &
                                        dJ123, dJ223, dJ323, &
                                        dJ133, dJ233, dJ333, -1_ik)

     call globalDiff2_gv (cctkGH, 0_ik, 1_ik, u, dxyu, J11, J21, J31, &
                                        J12, J22, J32, J13, J23, J33, &
                                        dJ111, dJ211, dJ311, &
                                        dJ112, dJ212, dJ312, &
                                        dJ113, dJ213, dJ313, &
                                        dJ122, dJ222, dJ322, &
                                        dJ123, dJ223, dJ323, &
                                        dJ133, dJ233, dJ333, -1_ik)

     call globalDiff2_gv (cctkGH, 0_ik, 2_ik, u, dxzu, J11, J21, J31, &
                                        J12, J22, J32, J13, J23, J33, &
                                        dJ111, dJ211, dJ311, &
                                        dJ112, dJ212, dJ312, &
                                        dJ113, dJ213, dJ313, &
                                        dJ122, dJ222, dJ322, &
                                        dJ123, dJ223, dJ323, &
                                        dJ133, dJ233, dJ333, -1_ik)

     call globalDiff2_gv (cctkGH, 1_ik, 1_ik, u, dyyu, J11, J21, J31, &
                                        J12, J22, J32, J13, J23, J33, &
                                        dJ111, dJ211, dJ311, &
                                        dJ112, dJ212, dJ312, &
                                        dJ113, dJ213, dJ313, &
                                        dJ122, dJ222, dJ322, &
                                        dJ123, dJ223, dJ323, &
                                        dJ133, dJ233, dJ333, -1_ik)

     call globalDiff2_gv (cctkGH, 1_ik, 2_ik, u, dyzu, J11, J21, J31, &
                                        J12, J22, J32, J13, J23, J33, &
                                        dJ111, dJ211, dJ311, &
                                        dJ112, dJ212, dJ312, &
                                        dJ113, dJ213, dJ313, &
                                        dJ122, dJ222, dJ322, &
                                        dJ123, dJ223, dJ323, &
                                        dJ133, dJ233, dJ333, -1_ik)

     call globalDiff2_gv (cctkGH, 2_ik, 2_ik, u, dzzu, J11, J21, J31, &
                                        J12, J22, J32, J13, J23, J33, &
                                        dJ111, dJ211, dJ311, &
                                        dJ112, dJ212, dJ312, &
                                        dJ113, dJ213, dJ313, &
                                        dJ122, dJ222, dJ322, &
                                        dJ123, dJ223, dJ323, &
                                        dJ133, dJ233, dJ333, -1_ik)
   end if

!!!$OMP PARALLEL WORKSHARE
   dxy_u = dxyu
   dxz_u = dxzu
   dyz_u = dyzu

   dxx_u = dxxu
   dyy_u = dyyu
   dzz_u = dzzu
   
   ! d/dt u = rho
   udot = rho
   
   ! d/dt rho = beta^i d_i rho + alpha / epsilon d_i vu^i
 !	write(*,*) "betax", betax
 !	write(*,*) "betay", betay
 !	write(*,*) "betaz", betaz
 !
 !	write (*,*) "dxwuxx", dxwuxx
 !	write (*,*) "epsilon", epsilon
 
 ! 	write (*,*) epsilon
 
   rhodot =   betax * dxrho + betay * dyrho + betaz * dzrho &
        &   + alpha / epsilon * &
              ( dxvux + dyvuy + dzvuz + &
                dxu * ( dxwuxx + dywuxy + dzwuxz) + &
                dyu * ( dxwuxy + dywuyy + dzwuyz) + &
                dzu * ( dxwuxz + dywuyz + dzwuzz) + &
                wuxx*dxxu + wuyy*dyyu + wuzz*dzzu + &
                2.0 * (wuxy*dxyu+wuxz*dxzu+wuyz*dyzu) )
 
!     do k=1,cctk_lsh(3)
!       do j=1, cctk_lsh(2)
!         do i=1, cctk_lsh(1)
!            rpt = sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
!            ! if (rpt<1.0) then
!              udot(i,j,k) = (1.0-1.0*exp(-rpt**2/0.7**2))*udot(i,j,k) - exp(-rpt**2/0.7**2)*u(i,j,k)
!              rhodot(i,j,k) = (1.0-1.0*exp(-rpt**2/0.7**2))*rhodot(i,j,k)  - exp(-rpt**2/0.7**2)*rho(i,j,k)
!            ! end if
!         end do
!       end do
!     end do
   
   ! d/dt v_i = d_i rho
   vxdot = dxrho
   vydot = dyrho
   vzdot = dzrho
!!!$OMP END PARALLEL WORKSHARE
 
   if (nonlinearrhs /= 0) then
!!!$OMP PARALLEL WORKSHARE
      rhodot = rhodot + u**powerrhs 
!!!$OMP END PARALLEL WORKSHARE
   end if
   
   deallocate ( dxrho, dyrho, dzrho, wuxx, wuxy, wuxz, wuyy, wuyz, wuzz, &
                dxwuxx, dxwuxy, dxwuxz, dxwuyy, dxwuyz, dxwuzz, &
                dywuxx, dywuxy, dywuxz, dywuyy, dywuyz, dywuzz, &
                dzwuxx, dzwuxy, dzwuxz, dzwuyy, dzwuyz, dzwuzz, &
                dxu, dyu, dzu, &
                dxxu, dxyu, dxzu, dyyu, dyzu, dzzu )

end subroutine LWT_calc_rhs
