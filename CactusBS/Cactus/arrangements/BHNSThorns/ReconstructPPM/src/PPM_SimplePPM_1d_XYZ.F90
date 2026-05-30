

 /*@@
   @file      PPM_SimplePPM_1d_XYZ.F90
   @date      Sun Feb 10 16:53:29 2002
   @author    Ian Hawke, Toni Font, Luca Baiotti, Frank Loeffler
   @desc 
   Routines to do PPM reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "PPM_Utils.h"

! let the preprocessor loop through x,y,z both with and without MHD

#if !defined(SimplePPM_1d_x_done)

#  define SimplePPM_1d_dir SimplePPM_1d_x
#  define BEGIN_LOOPS(o1,o2) do k = stencil,nrz-stencil+1 && do j = stencil,nry-stencil+1 && do i = (o1),nrx-(o2)
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i+(_i-i),j,k
#  define VIND3D(_i,_m) i+(_i-i),j,k,_m
#  define NTEMPS 16

#  define SimplePPM_1d_x_done

#elif !defined(SimplePPM_1d_y_done)

#  define SimplePPM_1d_dir SimplePPM_1d_y
#  define BEGIN_LOOPS(o1,o2) do k = stencil,nrz-stencil+1 && do j = (o1),nry-(o2) && do i = stencil,nrx-stencil+1
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i,j+(_i-i),k
#  define VIND3D(_i,_m) i,j+(_i-i),k,_m
#  define NTEMPS 16

#  define SimplePPM_1d_y_done

#elif !defined(SimplePPM_1d_z_done)

#  define SimplePPM_1d_dir SimplePPM_1d_z
#  define BEGIN_LOOPS(o1,o2) do k = (o1),nrz-(o2) && do j = stencil,nry-stencil+1 && do i = stencil,nrx-stencil+1
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i,j,k+(_i-i)
#  define VIND3D(_i,_m) i,j,k+(_i-i),_m
#  define NTEMPS 16

#  define SimplePPM_1d_z_done

#elif !defined(SimplePPM_1d_x_MHD_done)

#  define SimplePPM_1d_dir SimplePPM_1d_x_MHD
#  define BEGIN_LOOPS(o1,o2) do k = stencil,nrz-stencil+1 && do j = stencil,nry-stencil+1 && do i = (o1),nrx-(o2)
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i+(_i-i),j,k
#  define VIND3D(_i,_m) i+(_i-i),j,k,_m
#  define HAVE_BFIELD
#  define NTEMPS 22

#  define SimplePPM_1d_x_MHD_done

#elif !defined(SimplePPM_1d_y_MHD_done)

#  define SimplePPM_1d_dir SimplePPM_1d_y_MHD
#  define BEGIN_LOOPS(o1,o2) do k = stencil,nrz-stencil+1 && do j = (o1),nry-(o2) && do i = stencil,nrx-stencil+1
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i,j+(_i-i),k
#  define VIND3D(_i,_m) i,j+(_i-i),k,_m
#  define HAVE_BFIELD
#  define NTEMPS 22

#  define SimplePPM_1d_y_MHD_done

#elif !defined(SimplePPM_1d_z_MHD_done)

#  define SimplePPM_1d_dir SimplePPM_1d_z_MHD
#  define BEGIN_LOOPS(o1,o2) do k = (o1),nrz-(o2) && do j = stencil,nry-stencil+1 && do i = stencil,nrx-stencil+1
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i,j,k+(_i-i)
#  define VIND3D(_i,_m) i,j,k+(_i-i),_m
#  define HAVE_BFIELD
#  define NTEMPS 22

#  define SimplePPM_1d_z_MHD_done

#elif !defined(SimplePPM_1d_x_dMHD_done)

#  define SimplePPM_1d_dir SimplePPM_1d_x_dMHD
#  define BEGIN_LOOPS(o1,o2) do k = stencil,nrz-stencil+1 && do j = stencil,nry-stencil+1 && do i = (o1),nrx-(o2)
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i+(_i-i),j,k
#  define VIND3D(_i,_m) i+(_i-i),j,k,_m
#  define HAVE_BFIELD
#  define DIVCLEAN
#  define NTEMPS 24

#  define SimplePPM_1d_x_dMHD_done

#elif !defined(SimplePPM_1d_y_dMHD_done)

#  define SimplePPM_1d_dir SimplePPM_1d_y_dMHD
#  define BEGIN_LOOPS(o1,o2) do k = stencil,nrz-stencil+1 && do j = (o1),nry-(o2) && do i = stencil,nrx-stencil+1
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i,j+(_i-i),k
#  define VIND3D(_i,_m) i,j+(_i-i),k,_m
#  define HAVE_BFIELD
#  define DIVCLEAN
#  define NTEMPS 24

#  define SimplePPM_1d_y_dMHD_done

#elif !defined(SimplePPM_1d_z_dMHD_done)

#  define SimplePPM_1d_dir SimplePPM_1d_z_dMHD
#  define BEGIN_LOOPS(o1,o2) do k = (o1),nrz-(o2) && do j = stencil,nry-stencil+1 && do i = stencil,nrx-stencil+1
#  define END_LOOPS end do && end do && end do
#  define IND3D(_i) i,j,k+(_i-i)
#  define VIND3D(_i,_m) i,j,k+(_i-i),_m
#  define HAVE_BFIELD
#  define DIVCLEAN
#  define NTEMPS 24

#  define SimplePPM_1d_z_dMHD_done

#endif


! Start pre-processor loop here with module definition
#if defined(SimplePPM_1d_x_done) && !defined(SimplePPM_1d_y_done)
module PPM_SimplePPM_1d_XYZ
 CONTAINS
#endif


! Function arguments differ
#if !defined(HAVE_BFIELD)
  subroutine SimplePPM_1d_dir(handle,poly,&
     rho,velx,vely,velz,eps,press, &
     rhominus,velxminus,velyminus,velzminus,epsminus,&
     rhoplus,velxplus,velyplus,velzplus,epsplus,&
     trivial_rp,space_mask, excision_descriptors,&
     gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
     nrx, nry, nrz, ev_l, ev_r, xwind, stencil)
#elif !defined(DIVCLEAN)
  subroutine SimplePPM_1d_dir(handle,poly,&
     rho,velx,vely,velz,eps,press,bvec,bdirs, &
     rhominus,velxminus,velyminus,velzminus,epsminus,&
     bxminus,byminus,bzminus,&
     rhoplus,velxplus,velyplus,velzplus,epsplus,&
     bxplus,byplus,bzplus,&
     trivial_rp,space_mask, excision_descriptors,&
     gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
     nrx, nry, nrz, ev_l, ev_r, xwind, stencil)
#else
  subroutine SimplePPM_1d_dir(handle,poly,&
     rho,velx,vely,velz,eps,press,bvec,bdirs,dc_psi, &
     rhominus,velxminus,velyminus,velzminus,epsminus,&
     bxminus,byminus,bzminus,dc_psiminus,&
     rhoplus,velxplus,velyplus,velzplus,epsplus,&
     bxplus,byplus,bzplus,dc_psiplus,&
     trivial_rp,space_mask, excision_descriptors,&
     gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
     nrx, nry, nrz, ev_l, ev_r, xwind, stencil)
#endif

    USE PPM_Scalars
    USE PPM_EigenvalueCalcs

    implicit none

    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

    CCTK_INT, INTENT(IN) :: handle,poly,stencil
    CCTK_INT, INTENT(IN) :: nrx, nry, nrz
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(IN) :: rho,eps,press
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(IN), target :: velx,vely,velz
#ifdef HAVE_BFIELD 
    CCTK_REAL, dimension(nrx,nry,nrz,4), INTENT(IN) :: bvec
    CCTK_INT, dimension(3), INTENT(IN)            :: bdirs
# ifdef DIVCLEAN 
    CCTK_REAL, dimension(nrx,nry,nrz),  INTENT(IN) :: dc_psi
# endif
#endif

    CCTK_REAL, dimension(:,:,:), pointer :: vx,vy,vz
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(OUT) :: rhominus,velxminus,velyminus,velzminus,epsminus
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(OUT) :: rhoplus,velxplus,velyplus,velzplus,epsplus
#ifdef HAVE_BFIELD
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(OUT) :: bxminus,byminus,bzminus
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(OUT) :: bxplus,byplus,bzplus
# ifdef DIVCLEAN
    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(OUT) :: dc_psiminus,dc_psiplus
# endif
#endif


    CCTK_INT :: i,j,k,s
    CCTK_REAL, dimension(:,:,:), pointer :: drho,dvelx,dvely,dvelz,deps
    CCTK_REAL, dimension(:,:,:), pointer :: dmrho,dmvelx,dmvely,dmvelz,dmeps
    CCTK_REAL, dimension(:,:,:), pointer :: dpress,d2rho,tilde_flatten
#ifdef HAVE_BFIELD
    CCTK_REAL, dimension(:,:,:), pointer :: dbx,dby,dbz, dmbx,dmby,dmbz
# ifdef DIVCLEAN
    CCTK_REAL, dimension(:,:,:), pointer :: ddc_psi, dmdc_psi
# endif
#endif
    CCTK_REAL :: dpress2,dvel,w,flatten,eta,etatilde
    CCTK_REAL, dimension(:,:,:,:), target, allocatable :: temps

    logical, dimension(nrx,nry,nrz), INTENT(OUT) :: trivial_rp

    CCTK_INT, dimension(nrx,nry,nrz), INTENT(IN) :: space_mask
    CCTK_INT, dimension(3), INTENT(IN) :: excision_descriptors

    CCTK_REAL, dimension(nrx,nry,nrz), INTENT(IN) :: gxx, gxy, gxz, gyy, gyz, gzz, &
                                beta, alp, w_lorentz
    CCTK_REAL, dimension(nrx, nry, nrz), INTENT(OUT) :: ev_l, ev_r, xwind

    CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det
    CCTK_REAL, dimension(5) :: lam
    CCTK_REAL :: dupw, dloc, delta
    CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz
    integer :: ierr

    character(len=400) :: warnline

    !!$ get memory for the temporary variables
    allocate(temps(nrx,nry,nrz,NTEMPS), STAT=ierr)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Allocation problems with temps")
    end if
    !!$ get some nicer names for the temps
    vx => temps(1:nrx,1:nry,1:nrz,1)
    vy => temps(1:nrx,1:nry,1:nrz,2)
    vz => temps(1:nrx,1:nry,1:nrz,3)
    drho => temps(1:nrx,1:nry,1:nrz,4)
    dvelx => temps(1:nrx,1:nry,1:nrz,5)
    dvely => temps(1:nrx,1:nry,1:nrz,6)
    dvelz => temps(1:nrx,1:nry,1:nrz,7)
    deps => temps(1:nrx,1:nry,1:nrz,8)
    dmrho => temps(1:nrx,1:nry,1:nrz,9)
    dmvelx => temps(1:nrx,1:nry,1:nrz,10)
    dmvely => temps(1:nrx,1:nry,1:nrz,11)
    dmvelz => temps(1:nrx,1:nry,1:nrz,12)
    dmeps => temps(1:nrx,1:nry,1:nrz,13)
    dpress => temps(1:nrx,1:nry,1:nrz,14)
    d2rho => temps(1:nrx,1:nry,1:nrz,15)
    tilde_flatten => temps(1:nrx,1:nry,1:nrz,16)
#ifdef HAVE_BFIELD
    dbx => temps(1:nrx,1:nry,1:nrz,17)
    dby => temps(1:nrx,1:nry,1:nrz,18)
    dbz => temps(1:nrx,1:nry,1:nrz,19)
    dmbx => temps(1:nrx,1:nry,1:nrz,20)
    dmby => temps(1:nrx,1:nry,1:nrz,21)
    dmbz => temps(1:nrx,1:nry,1:nrz,22)
# ifdef DIVCLEAN
    ddc_psi => temps(1:nrx,1:nry,1:nrz,23)
    dmdc_psi => temps(1:nrx,1:nry,1:nrz,24)
# endif
#endif

    !$OMP PARALLEL &
    !$OMP PRIVATE (i,j,k,s,dpress2,dvel,w,flatten,eta,etatilde, &
    !$OMP          uxx,uxy,uxz,uyy,uyz,uzz, det, agxx,agxy,agxz,agyy,agyz,agzz, &
    !$OMP          lam, dupw, dloc, delta, warnline)

!!$  Initially, all the Riemann problems will be trivial
    !$OMP WORKSHARE
    trivial_rp = .true.
    !$OMP END WORKSHARE NOWAIT

!!$ reconstruct w^i = w_lorentz*v^i to ensure slower than light speeds?
    if (reconstruct_Wv.ne.0) then
      !$OMP WORKSHARE
      vx = w_lorentz*velx
      vy = w_lorentz*vely
      vz = w_lorentz*velz
      !$OMP END WORKSHARE NOWAIT
    else
      !$OMP SINGLE
      vx => velx
      vy => vely
      vz => velz
      !$OMP END SINGLE NOWAIT
    end if

    if(ppm_mppm.ne.0) then

      !$OMP WORKSHARE 
      ev_l=0.d0
      ev_r=0.d0
      xwind=0.d0
      !$OMP END WORKSHARE NOWAIT

    end if
    !$OMP BARRIER

!!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
!!$  This is the expression for an even grid.

    !$OMP DO
    BEGIN_LOOPS(2,1)
      drho(IND3D(i)) = 0.5d0 * (rho(IND3D(i+1)) - rho(IND3D(i-1)))
      d2rho(IND3D(i)) = (rho(IND3D(i+1)) - 2.d0 * rho(IND3D(i)) + rho(IND3D(i-1)))! / 6.d0 / dx / dx
      ! since we use d2rho only for the condition d2rho(IND3D(i+1))*d2rhoi(IND3D(i-1))<0 
      ! the denominator is not necessary
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,1)
      dvelx(IND3D(i)) = 0.5d0 * (vx(IND3D(i+1)) - vx(IND3D(i-1)))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,1)
      dvely(IND3D(i)) = 0.5d0 * (vy(IND3D(i+1)) - vy(IND3D(i-1)))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,1)
      dvelz(IND3D(i)) = 0.5d0 * (vz(IND3D(i+1)) - vz(IND3D(i-1)))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,1)
      dpress(IND3D(i)) = press(IND3D(i+1)) - press(IND3D(i-1))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,1)
      deps(IND3D(i)) = 0.5d0 * (eps(IND3D(i+1)) - eps(IND3D(i-1)))
    END_LOOPS
    !$OMP END DO NOWAIT
#if defined(HAVE_BFIELD)
# if defined(DIVCLEAN)
    !$OMP DO
    BEGIN_LOOPS(2,1)
      dbx(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i+1,bdirs(1))) - bvec(VIND3D(i-1,bdirs(1))))
      dby(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i+1,bdirs(2))) - bvec(VIND3D(i-1,bdirs(2))))
      dbz(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i+1,bdirs(3))) - bvec(VIND3D(i-1,bdirs(3))))
      ddc_psi(IND3D(i)) = 0.5d0 * (dc_psi(IND3D(i+1)) - dc_psi(IND3D(i-1)))
    END_LOOPS
    !$OMP END DO NOWAIT
# else
    !$OMP DO
    BEGIN_LOOPS(2,1)
      dbx(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i+1,bdirs(1))) - bvec(VIND3D(i-1,bdirs(1))))
      dby(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i+1,bdirs(2))) - bvec(VIND3D(i-1,bdirs(2))))
      dbz(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i+1,bdirs(3))) - bvec(VIND3D(i-1,bdirs(3))))
    END_LOOPS
    !$OMP END DO NOWAIT
# endif
#endif

    !$OMP BARRIER

!!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

#define STEEP(x,dx,dmx)                                              \
    BEGIN_LOOPS(2,1)                                                 &&\
      if ( (x(IND3D(i+1)) - x(IND3D(i))) * (x(IND3D(i)) - x(IND3D(i-1))) > 0.d0 ) then           &&\
        dmx(IND3D(i)) = sign(1.d0, dx(IND3D(i))) *                           \
             min(abs(dx(IND3D(i))), 2.d0 * abs(x(IND3D(i)) - x(IND3D(i-1))),    \
                             2.d0 * abs(x(IND3D(i+1)) - x(IND3D(i))))      &&\
      else                                                           &&\
        dmx(IND3D(i)) = 0.d0                                            &&\
      end if                                                         &&\
    END_LOOPS
#define VSTEEP(x,dx,dmx,bi)                                              \
    BEGIN_LOOPS(2,1)                                                 &&\
      if ( (x(VIND3D(i+1,bi)) - x(VIND3D(i,bi))) * (x(VIND3D(i,bi)) - x(VIND3D(i-1,bi))) > 0.d0 ) then           &&\
        dmx(IND3D(i)) = sign(1.d0, dx(IND3D(i))) *                           \
             min(abs(dx(IND3D(i))), 2.d0 * abs(x(VIND3D(i,bi)) - x(VIND3D(i-1,bi))),    \
                             2.d0 * abs(x(VIND3D(i+1,bi)) - x(VIND3D(i,bi))))      &&\
      else                                                           &&\
        dmx(IND3D(i)) = 0.d0                                            &&\
      end if                                                         &&\
    END_LOOPS
    !$OMP DO SCHEDULE(dynamic)
    STEEP(rho, drho, dmrho)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    STEEP(vx, dvelx, dmvelx)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    STEEP(vy, dvely, dmvely)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    STEEP(vz, dvelz, dmvelz)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    STEEP(eps, deps, dmeps)
    !$OMP END DO NOWAIT
#if defined(HAVE_BFIELD)
    !$OMP DO SCHEDULE(dynamic)
    VSTEEP(bvec, dbx, dmbx, bdirs(1))
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    VSTEEP(bvec, dby, dmby, bdirs(2))
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    VSTEEP(bvec, dbz, dmbz, bdirs(3))
    !$OMP END DO NOWAIT
# if defined(DIVCLEAN)
    !$OMP DO SCHEDULE(dynamic)
    STEEP(dc_psi, ddc_psi, dmdc_psi)
    !$OMP END DO NOWAIT
# endif
#endif
    !$OMP BARRIER

!!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

    !$OMP DO
    BEGIN_LOOPS(2,2)
      rhoplus(IND3D(i)) = 0.5d0 * (rho(IND3D(i)) + rho(IND3D(i+1))) + &
           (dmrho(IND3D(i)) - dmrho(IND3D(i+1))) / 6.d0
      rhominus(IND3D(i+1)) = rhoplus(IND3D(i))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,2)
      velxplus(IND3D(i)) = 0.5d0 * (vx(IND3D(i)) + vx(IND3D(i+1))) + &
           (dmvelx(IND3D(i)) - dmvelx(IND3D(i+1))) / 6.d0
      velxminus(IND3D(i+1)) = velxplus(IND3D(i))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,2)
      velyplus(IND3D(i)) = 0.5d0 * (vy(IND3D(i)) + vy(IND3D(i+1))) + &
           (dmvely(IND3D(i)) - dmvely(IND3D(i+1))) / 6.d0
      velyminus(IND3D(i+1)) = velyplus(IND3D(i))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,2)
      velzplus(IND3D(i)) = 0.5d0 * (vz(IND3D(i)) + vz(IND3D(i+1))) + &
           (dmvelz(IND3D(i)) - dmvelz(IND3D(i+1))) / 6.d0
      velzminus(IND3D(i+1)) = velzplus(IND3D(i))
    END_LOOPS
    !$OMP END DO NOWAIT
    !$OMP DO
    BEGIN_LOOPS(2,2)
      epsplus(IND3D(i)) = 0.5d0 * (eps(IND3D(i)) + eps(IND3D(i+1))) + &
           (dmeps(IND3D(i)) - dmeps(IND3D(i+1))) / 6.d0
      epsminus(IND3D(i+1)) = epsplus(IND3D(i))
    END_LOOPS
    !$OMP END DO NOWAIT
#if defined(HAVE_BFIELD)
    !$OMP DO
    BEGIN_LOOPS(2,2)
      bxplus(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i,bdirs(1))) + bvec(VIND3D(i+1,bdirs(1)))) + &
            (dmbx(IND3D(i)) - dmbx(IND3D(i+1))) / 6.d0
      bxminus(IND3D(i+1)) = bxplus(IND3D(i))
      byplus(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i,bdirs(2))) + bvec(VIND3D(i+1,bdirs(2)))) + &
            (dmby(IND3D(i)) - dmby(IND3D(i+1))) / 6.d0
      byminus(IND3D(i+1)) = byplus(IND3D(i))
      bzplus(IND3D(i)) = 0.5d0 * (bvec(VIND3D(i,bdirs(3))) + bvec(VIND3D(i+1,bdirs(3)))) + &
            (dmbz(IND3D(i)) - dmbz(IND3D(i+1))) / 6.d0
      bzminus(IND3D(i+1)) = bzplus(IND3D(i))
# if defined(DIVCLEAN)
      dc_psiplus(IND3D(i)) = 0.5d0 * (dc_psi(IND3D(i)) + dc_psi(IND3D(i+1))) + &
           (dmdc_psi(IND3D(i)) - dmdc_psi(IND3D(i+1))) / 6.d0
      dc_psiminus(IND3D(i+1)) = dc_psiplus(IND3D(i))
# endif
    END_LOOPS
    !$OMP END DO NOWAIT
#endif
    !$OMP BARRIER


!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which may be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even. 
!!$Note also that we don''t have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).

    if (ppm_detect .ne. 0) then

      !$OMP DO SCHEDULE(dynamic)
      BEGIN_LOOPS(3,2)
        if ( (d2rho(IND3D(i+1))*d2rho(IND3D(i-1)) < 0.d0).and.(abs(rho(IND3D(i+1))-rho(IND3D(i-1))) - &
             ppm_epsilon_shock * min(abs(rho(IND3D(i+1))), abs(rho(IND3D(i-1)))) > 0.d0) ) then
          etatilde = (rho(IND3D(i-2)) - rho(IND3D(i+2)) + 4.d0 * drho(IND3D(i))) / (drho(IND3D(i)) * 12.d0)
        else
          etatilde = 0.d0
        end if
        eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
        if (ppm_k0 * abs(drho(IND3D(i))) * min(press(IND3D(i-1)),press(IND3D(i+1))) < &
             abs(dpress(IND3D(i))) * min(rho(IND3D(i-1)), rho(IND3D(i+1)))) then
          eta = 0.d0
        end if
        if (eta > 0.d0) then
          trivial_rp(IND3D(i-1)) = .false.
          trivial_rp(IND3D(i)) = .false.
        end if
        rhominus(IND3D(i)) = rhominus(IND3D(i)) * (1.d0 - eta) + &
             (rho(IND3D(i-1)) + 0.5d0 * dmrho(IND3D(i-1))) * eta
        rhoplus(IND3D(i)) = rhoplus(IND3D(i)) * (1.d0 - eta) + &
             (rho(IND3D(i+1)) - 0.5d0 * dmrho(IND3D(i+1))) * eta
      END_LOOPS
      !$OMP END DO 

    end if

    !!$ mppm
#define D_UPW(x) (0.5d0 * (x(IND3D(i)) + x(IND3D(i+1))))
#define LEFT1(x)  (13.d0*x(IND3D(i+1))-5.d0*x(IND3D(i+2))+x(IND3D(i+3))+3.d0*x(IND3D(i  )))/12.d0
#define RIGHT1(x) (13.d0*x(IND3D(i  ))-5.d0*x(IND3D(i-1))+x(IND3D(i-2))+3.d0*x(IND3D(i+1)))/12.d0
#define VD_UPW(x,bi) (0.5d0 * (x(VIND3D(i,bi)) + x(VIND3D(i+1,bi))))
#define VLEFT1(x,bi) (13.d0*x(VIND3D(i+1,bi))-5.d0*x(VIND3D(i+2,bi))+x(VIND3D(i+3,bi))+3.d0*x(VIND3D(i  ,bi)))/12.d0
#define VRIGHT1(x,bi)(13.d0*x(VIND3D(i  ,bi))-5.d0*x(VIND3D(i-1,bi))+x(VIND3D(i-2,bi))+3.d0*x(VIND3D(i+1,bi)))/12.d0
    if (ppm_mppm .gt. 0) then
      !$OMP DO SCHEDULE(dynamic)
      BEGIN_LOOPS(3,3)
        agxx = 0.5d0*( gxx(IND3D(i)) + gxx(IND3D(i+1)) )
        agxy = 0.5d0*( gxy(IND3D(i)) + gxy(IND3D(i+1)) )
        agxz = 0.5d0*( gxz(IND3D(i)) + gxz(IND3D(i+1)) )
        agyy = 0.5d0*( gyy(IND3D(i)) + gyy(IND3D(i+1)) )
        agyz = 0.5d0*( gyz(IND3D(i)) + gyz(IND3D(i+1)) )
        agzz = 0.5d0*( gzz(IND3D(i)) + gzz(IND3D(i+1)) )
        det=SPATIAL_DET(agxx,agxy,agxz,agyy,agyz,agzz)
        call UpperMetric (uxx, uxy, uxz, uyy, uyz, uzz, &
                          det, agxx, agxy, agxz, agyy, agyz, agzz)
!!$ these always have to be vel[xyz], not v[xyz] since they are used to
!!$ calculate the eigenvalues
        call eigenvalues_ppm(handle,&
                   D_UPW(rho), D_UPW(velx), D_UPW(vely), D_UPW(velz), &
                   D_UPW(eps), D_UPW(w_lorentz), lam, &
                   agxx, agxy, agxz, agyy, agyz, agzz, &
                   uxx, det, D_UPW(alp), D_UPW(beta))
        ev_l(IND3D(i))=lam(1)
        ev_r(IND3D(i))=lam(5)
        xwind(IND3D(i)) = (lam(1) + lam(5)) / (abs(lam(1)) + abs(lam(5)))
        xwind(IND3D(i)) = min(1.d0, max(-1.d0, xwind(IND3D(i))))
#define LEFTPLUS(x,xplus)    xplus(IND3D(i))   =       abs(xwind(IND3D(i)))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xplus(IND3D(i))
#define LEFTMINUS(x,xminus)  xminus(IND3D(i+1))=       abs(xwind(IND3D(i)))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xminus(IND3D(i+1))
#define RIGHTPLUS(x,xplus)   xplus(IND3D(i))   =       abs(xwind(IND3D(i)))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xplus(IND3D(i))
#define RIGHTMINUS(x,xminus) xminus(IND3D(i+1))=       abs(xwind(IND3D(i)))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xminus(IND3D(i+1))
#define CHECK(x,xc) if (x(IND3D(i+1)) .gt. x(IND3D(i))) then && xc=min(x(IND3D(i+1)),max(x(IND3D(i)),xc)) && else && xc=min(x(IND3D(i)),max(x(IND3D(i+1)),xc)) && endif
#define VLEFTPLUS(x,xplus,bi)    xplus(IND3D(i))   =       abs(xwind(IND3D(i)))  * VLEFT1(x,bi) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xplus(IND3D(i))
#define VLEFTMINUS(x,xminus,bi)  xminus(IND3D(i+1))=       abs(xwind(IND3D(i)))  * VLEFT1(x,bi) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xminus(IND3D(i+1))
#define VRIGHTPLUS(x,xplus,bi)   xplus(IND3D(i))   =       abs(xwind(IND3D(i)))  * VRIGHT1(x,bi) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xplus(IND3D(i))
#define VRIGHTMINUS(x,xminus,bi) xminus(IND3D(i+1))=       abs(xwind(IND3D(i)))  * VRIGHT1(x,bi) + \
                                          (1.d0-abs(xwind(IND3D(i)))) * xminus(IND3D(i+1))
#define VCHECK(x,xc,bi) if (x(VIND3D(i+1,bi)) .gt. x(VIND3D(i,bi))) then && xc=min(x(VIND3D(i+1,bi)),max(x(VIND3D(i,bi)),xc)) && else && xc=min(x(VIND3D(i,bi)),max(x(VIND3D(i+1,bi)),xc)) && endif
!!$      xwind(i)=0.d0
        if (xwind(IND3D(i)) .lt. 0.0d0) then
          LEFTPLUS(rho, rhoplus)
          LEFTMINUS(rho, rhominus)
          LEFTPLUS(vx, velxplus)
          LEFTMINUS(vx, velxminus)
          LEFTPLUS(vy, velyplus)
          LEFTMINUS(vy, velyminus)
          LEFTPLUS(vz, velzplus)
          LEFTMINUS(vz, velzminus)
          LEFTPLUS(eps, epsplus)
          LEFTMINUS(eps, epsminus)
#if defined(HAVE_BFIELD)
          VLEFTPLUS(bvec, bxplus, bdirs(1))
          VLEFTMINUS(bvec, bxminus, bdirs(1))
          VLEFTPLUS(bvec, byplus, bdirs(2))
          VLEFTMINUS(bvec, byminus, bdirs(2))
          VLEFTPLUS(bvec, bzplus, bdirs(3))
          VLEFTMINUS(bvec, bzminus, bdirs(3))
# if defined(DIVCLEAN)
          LEFTPLUS(dc_psi, dc_psiplus)
          LEFTPLUS(dc_psi, dc_psiminus)
# endif
#endif
        else
          RIGHTPLUS(rho, rhoplus)
          RIGHTMINUS(rho, rhominus)
          RIGHTPLUS(vx, velxplus)
          RIGHTMINUS(vx, velxminus)
          RIGHTPLUS(vy, velyplus)
          RIGHTMINUS(vy, velyminus)
          RIGHTPLUS(vz, velzplus)
          RIGHTMINUS(vz, velzminus)
          RIGHTPLUS(eps, epsplus)
          RIGHTMINUS(eps, epsminus)
#if defined(HAVE_BFIELD)
          VRIGHTPLUS(bvec, bxplus, bdirs(1))
          VRIGHTMINUS(bvec, bxminus, bdirs(1))
          VRIGHTPLUS(bvec, byplus, bdirs(2))
          VRIGHTMINUS(bvec, byminus, bdirs(2))
          VRIGHTPLUS(bvec, bzplus, bdirs(3))
          VRIGHTMINUS(bvec, bzminus, bdirs(3))
# if defined(DIVCLEAN)
          RIGHTPLUS(dc_psi, dc_psiplus)
          RIGHTPLUS(dc_psi, dc_psiminus)
# endif
#endif
        end if
        CHECK(rho, rhoplus(IND3D(i)))
        CHECK(rho, rhominus(IND3D(i+1)))
        CHECK(vx, velxplus(IND3D(i)))
        CHECK(vx, velxminus(IND3D(i+1)))
        CHECK(vy, velyplus(IND3D(i)))
        CHECK(vy, velyminus(IND3D(i+1)))
        CHECK(vz, velzplus(IND3D(i)))
        CHECK(vz, velzminus(IND3D(i+1)))
        CHECK(eps, epsplus(IND3D(i)))
        CHECK(eps, epsminus(IND3D(i+1)))
#if defined(HAVE_BFIELD)
        VCHECK(bvec, bxplus(IND3D(i)), bdirs(1))
        VCHECK(bvec, bxminus(IND3D(i+1)), bdirs(1))
        VCHECK(bvec, byplus(IND3D(i)), bdirs(2))
        VCHECK(bvec, byminus(IND3D(i+1)), bdirs(2))
        VCHECK(bvec, bzplus(IND3D(i)), bdirs(3))
        VCHECK(bvec, bzminus(IND3D(i+1)), bdirs(3))
# if defined(DIVCLEAN)
        CHECK(dc_psi, dc_psiplus(IND3D(i)))
        CHECK(dc_psi, dc_psiminus(IND3D(i+1)))
# endif
#endif
!!$      if ((j .eq. 4) .and. (k .eq. 4)) then
!!$        write (*,*) rhoplus(IND3D(i)), rhominus(IND3D(i+1))
!!$      end if
      END_LOOPS
      !$OMP END DO
    end if

!!$  Zone flattening. See appendix of C&W, p. 197-8.

    !$OMP DO SCHEDULE(dynamic)
    BEGIN_LOOPS(3,2)
      dpress2 = press(IND3D(i+2)) - press(IND3D(i-2))
      dvel = vx(IND3D(i+1)) - vx(IND3D(i-1))
      if ( (abs(dpress(IND3D(i))) >  ppm_epsilon * min(press(IND3D(i-1)),press(IND3D(i+1)))) .and. &
           (dvel < 0.d0) .and. abs(dpress2) >= ppm_small) then
        tilde_flatten(IND3D(i)) = max(0.d0, 1.d0 - max(0.d0, ppm_omega2 * &
             (dpress(IND3D(i)) / dpress2 - ppm_omega1)))
      else
        tilde_flatten(IND3D(i)) = 1.d0
      end if
      if (tilde_flatten(IND3D(i)) .ne. 1.d0) then
        trivial_rp(IND3D(i-1)) = .false.
        trivial_rp(IND3D(i)) = .false.
      end if
    END_LOOPS
    !$OMP END DO




    if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
      !$OMP DO
      BEGIN_LOOPS(3,2)
        flatten = tilde_flatten(IND3D(i))
        rhoplus(IND3D(i)) = flatten * rhoplus(IND3D(i)) + (1.d0 - flatten) * rho(IND3D(i))
        rhominus(IND3D(i)) = flatten * rhominus(IND3D(i)) + (1.d0 - flatten) * rho(IND3D(i))

        velxplus(IND3D(i)) = flatten * velxplus(IND3D(i)) + (1.d0 - flatten) * vx(IND3D(i))
        velxminus(IND3D(i)) = flatten * velxminus(IND3D(i)) + (1.d0 - flatten) * vx(IND3D(i))

        velyplus(IND3D(i)) = flatten * velyplus(IND3D(i)) + (1.d0 - flatten) * vy(IND3D(i))
        velyminus(IND3D(i)) = flatten * velyminus(IND3D(i)) + (1.d0 - flatten) * vy(IND3D(i))
 
        velzplus(IND3D(i)) = flatten * velzplus(IND3D(i)) + (1.d0 - flatten) * vz(IND3D(i))
        velzminus(IND3D(i)) = flatten * velzminus(IND3D(i)) + (1.d0 - flatten) * vz(IND3D(i))

        epsplus(IND3D(i)) = flatten * epsplus(IND3D(i)) + (1.d0 - flatten) * eps(IND3D(i))
        epsminus(IND3D(i)) = flatten * epsminus(IND3D(i)) + (1.d0 - flatten) * eps(IND3D(i))

#if defined(HAVE_BFIELD)
        bxplus(IND3D(i)) = flatten * bxplus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(1)))
        bxminus(IND3D(i)) = flatten * bxminus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(1)))
        byplus(IND3D(i)) = flatten * byplus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(2)))
        byminus(IND3D(i)) = flatten * byminus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(2)))
        bzplus(IND3D(i)) = flatten * bzplus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(3)))
        bzminus(IND3D(i)) = flatten * bzminus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(3)))
# if defined(DIVCLEAN)
        dc_psiplus(IND3D(i)) = flatten * dc_psiplus(IND3D(i)) + (1.d0 - flatten) * dc_psi(IND3D(i))
        dc_psiminus(IND3D(i)) = flatten * dc_psiminus(IND3D(i)) + (1.d0 - flatten) * dc_psi(IND3D(i))
# endif
#endif

      END_LOOPS
      !$OMP END DO

    else  !!$ Really implement C&W, page 197; which requires stencil 4.
      !$OMP DO
      BEGIN_LOOPS(4,3)
        s=sign(1.d0, -dpress(IND3D(i)))
        flatten = max(tilde_flatten(IND3D(i)), tilde_flatten(IND3D(i+s)))  
        rhoplus(IND3D(i)) = flatten * rhoplus(IND3D(i)) + (1.d0 - flatten) * rho(IND3D(i))
        rhominus(IND3D(i)) = flatten * rhominus(IND3D(i)) + (1.d0 - flatten) * rho(IND3D(i))
        velxplus(IND3D(i)) = flatten * velxplus(IND3D(i)) + (1.d0 - flatten) * vx(IND3D(i))
        velxminus(IND3D(i)) = flatten * velxminus(IND3D(i)) + (1.d0 - flatten) * vx(IND3D(i))
        velyplus(IND3D(i)) = flatten * velyplus(IND3D(i)) + (1.d0 - flatten) * vy(IND3D(i))
        velyminus(IND3D(i)) = flatten * velyminus(IND3D(i)) + (1.d0 - flatten) * vy(IND3D(i))
        velzplus(IND3D(i)) = flatten * velzplus(IND3D(i)) + (1.d0 - flatten) * vz(IND3D(i))
        velzminus(IND3D(i)) = flatten * velzminus(IND3D(i)) + (1.d0 - flatten) * vz(IND3D(i))
        epsplus(IND3D(i)) = flatten * epsplus(IND3D(i)) + (1.d0 - flatten) * eps(IND3D(i))
        epsminus(IND3D(i)) = flatten * epsminus(IND3D(i)) + (1.d0 - flatten) * eps(IND3D(i))
#if defined(HAVE_BFIELD)
        bxplus(IND3D(i)) = flatten * bxplus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(1)))
        bxminus(IND3D(i)) = flatten * bxminus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(1)))
        byplus(IND3D(i)) = flatten * byplus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(2)))
        byminus(IND3D(i)) = flatten * byminus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(2)))
        bzplus(IND3D(i)) = flatten * bzplus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(3)))
        bzminus(IND3D(i)) = flatten * bzminus(IND3D(i)) + (1.d0 - flatten) * bvec(VIND3D(i,bdirs(3)))
# if defined(DIVCLEAN)
        dc_psiplus(IND3D(i)) = flatten * dc_psiplus(IND3D(i)) + (1.d0 - flatten) * dc_psi(IND3D(i))
        dc_psiminus(IND3D(i)) = flatten * dc_psiminus(IND3D(i)) + (1.d0 - flatten) * dc_psi(IND3D(i))
# endif
#endif
      END_LOOPS
      !$OMP END DO
    end if


!!$ Monotonicity. See (1.10) of C&W.

#define MON(xminus,x,xplus)                                           \
    BEGIN_LOOPS(stencil, stencil - 1)                                 &&\
      if (.not.( (xplus(IND3D(i)).eq.x(IND3D(i))) .and. (x(IND3D(i)).eq.xminus(IND3D(i))) )     \
          .and. ((xplus(IND3D(i))-x(IND3D(i)))*(x(IND3D(i))-xminus(IND3D(i))) .le. 0.d0)) then&&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
        xminus(IND3D(i)) = x(IND3D(i))                                      &&\
        xplus(IND3D(i)) = x(IND3D(i))                                       &&\
      else if (6.d0 * (xplus(IND3D(i)) - xminus(IND3D(i))) * (x(IND3D(i)) - 0.5d0 *      \
                  (xplus(IND3D(i)) + xminus(IND3D(i)))) >                     \
                  (xplus(IND3D(i)) - xminus(IND3D(i)))**2) then             &&\
        xminus(IND3D(i)) = 3.d0 * x(IND3D(i)) - 2.d0 * xplus(IND3D(i))         &&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
      else if (6.d0 * (xplus(IND3D(i)) - xminus(IND3D(i))) * (x(IND3D(i)) - 0.5d0 *      \
                  (xplus(IND3D(i)) + xminus(IND3D(i)))) <                     \
                 -(xplus(IND3D(i)) - xminus(IND3D(i)))**2) then             &&\
        xplus(IND3D(i)) = 3.d0 * x(IND3D(i)) - 2.d0 * xminus(IND3D(i))         &&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
      end if                                                          &&\
      if (.not.( (xplus(IND3D(i)).eq.x(IND3D(i))) .and. (x(IND3D(i)).eq.xminus(IND3D(i))) ) ) then     &&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
      end if                                                          &&\
    END_LOOPS
#define VMON(xminus,x,xplus,bi)                                           \
    BEGIN_LOOPS(stencil, stencil - 1)                                 &&\
      if (.not.( (xplus(IND3D(i)).eq.x(VIND3D(i,bi))) .and. (x(VIND3D(i,bi)).eq.xminus(IND3D(i))) )     \
          .and. ((xplus(IND3D(i))-x(VIND3D(i,bi)))*(x(VIND3D(i,bi))-xminus(IND3D(i))) .le. 0.d0)) then&&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
        xminus(IND3D(i)) = x(VIND3D(i,bi))                                      &&\
        xplus(IND3D(i)) = x(VIND3D(i,bi))                                       &&\
      else if (6.d0 * (xplus(IND3D(i)) - xminus(IND3D(i))) * (x(VIND3D(i,bi)) - 0.5d0 *      \
                  (xplus(IND3D(i)) + xminus(IND3D(i)))) >                     \
                  (xplus(IND3D(i)) - xminus(IND3D(i)))**2) then             &&\
        xminus(IND3D(i)) = 3.d0 * x(VIND3D(i,bi)) - 2.d0 * xplus(IND3D(i))         &&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
      else if (6.d0 * (xplus(IND3D(i)) - xminus(IND3D(i))) * (x(VIND3D(i,bi)) - 0.5d0 *      \
                  (xplus(IND3D(i)) + xminus(IND3D(i)))) <                     \
                 -(xplus(IND3D(i)) - xminus(IND3D(i)))**2) then             &&\
        xplus(IND3D(i)) = 3.d0 * x(VIND3D(i,bi)) - 2.d0 * xminus(IND3D(i))         &&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
      end if                                                          &&\
      if (.not.( (xplus(IND3D(i)).eq.x(VIND3D(i,bi))) .and. (x(VIND3D(i,bi)).eq.xminus(IND3D(i))) ) ) then     &&\
        trivial_rp(IND3D(i-1)) = .false.                                 &&\
        trivial_rp(IND3D(i)) = .false.                                   &&\
      end if                                                          &&\
    END_LOOPS



    !$OMP DO SCHEDULE(dynamic)
    MON(rhominus,rho,rhoplus)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    MON(velxminus,vx,velxplus)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    MON(velyminus,vy,velyplus)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    MON(velzminus,vz,velzplus)
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    MON(epsminus,eps,epsplus)
    !$OMP END DO NOWAIT
#if defined(HAVE_BFIELD)
    !$OMP DO SCHEDULE(dynamic)
    VMON(bxminus,bvec,bxplus,bdirs(1))
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    VMON(byminus,bvec,byplus,bdirs(2))
    !$OMP END DO NOWAIT
    !$OMP DO SCHEDULE(dynamic)
    VMON(bzminus,bvec,bzplus,bdirs(3))
    !$OMP END DO NOWAIT
# if defined(DIVCLEAN)
    !$OMP DO SCHEDULE(dynamic)
    MON(dc_psiminus,dc_psi,dc_psiplus)
    !$OMP END DO NOWAIT
# endif
#endif
    !$OMP BARRIER

    !!$ excision

    !!$ recover primitive v^i from W*v^i at interfaces
    if (reconstruct_Wv.ne.0) then
      !$OMP DO
      BEGIN_LOOPS(stencil, stencil - 1)
        agxx = 0.5d0*( gxx(IND3D(i)) + gxx(IND3D(i+1)) )
        agxy = 0.5d0*( gxy(IND3D(i)) + gxy(IND3D(i+1)) )
        agxz = 0.5d0*( gxz(IND3D(i)) + gxz(IND3D(i+1)) )
        agyy = 0.5d0*( gyy(IND3D(i)) + gyy(IND3D(i+1)) )
        agyz = 0.5d0*( gyz(IND3D(i)) + gyz(IND3D(i+1)) )
        agzz = 0.5d0*( gzz(IND3D(i)) + gzz(IND3D(i+1)) )
        ! we re-use the variable w defined for zone flattening to hold the Lorentz
        ! factor
        ! divide out the Lorentz factor for both the plus and minus quantities
        ! this should by construction ensure that any Lorentz factor calculated
        ! from them later on is physical (ie. > 1.d0)
        w = sqrt( 1.d0 + agxx*velxminus(IND3D(i))*velxminus(IND3D(i)) + agyy*velyminus(IND3D(i))*velyminus(IND3D(i)) &
                + agzz*velzminus(IND3D(i))*velzminus(IND3D(i)) + 2.d0*agxy*velxminus(IND3D(i))*velyminus(IND3D(i)) &
                + 2.d0*agxz*velxminus(IND3D(i))*velzminus(IND3D(i)) + 2.d0*agyz*velyminus(IND3D(i))*velzminus(IND3D(i)) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(IND3D(i)).eq.1.d0 (no CCTK_ARGUMENTS)
          !$OMP CRITICAL
          write(warnline,'("SimplePPM_1d (minus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxminus(IND3D(i)),velyminus(IND3D(i)),velzminus(IND3D(i))
          call CCTK_WARN(1,warnline)
          !$OMP END CRITICAL
        end if 
        velxminus(IND3D(i)) = velxminus(IND3D(i))/w
        velyminus(IND3D(i)) = velyminus(IND3D(i))/w
        velzminus(IND3D(i)) = velzminus(IND3D(i))/w
        w = sqrt( 1.d0 + agxx*velxplus(IND3D(i))*velxplus(IND3D(i)) + agyy*velyplus(IND3D(i))*velyplus(IND3D(i)) &
                + agzz*velzplus(IND3D(i))*velzplus(IND3D(i)) + 2.d0*agxy*velxplus(IND3D(i))*velyplus(IND3D(i)) &
                + 2.d0*agxz*velxplus(IND3D(i))*velzplus(IND3D(i)) + 2.d0*agyz*velyplus(IND3D(i))*velzplus(IND3D(i)) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(IND3D(i)).eq.1.d0 (no CCTK_ARGUMENTS)
          !$OMP CRITICAL
          write(warnline,'("SimplePPM_1d (plus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxplus(IND3D(i)),velyplus(IND3D(i)),velzplus(IND3D(i))
          call CCTK_WARN(1,warnline)
          !$OMP END CRITICAL
        end if 
        velxplus(IND3D(i)) = velxplus(IND3D(i))/w
        velyplus(IND3D(i)) = velyplus(IND3D(i))/w
        velzplus(IND3D(i)) = velzplus(IND3D(i))/w
      END_LOOPS
      !$OMP END DO
    end if

    if (check_for_trivial_rp .eq. 0) then
      !$OMP WORKSHARE
      trivial_rp = .false.
      !$OMP END WORKSHARE
    end if
    !$OMP END PARALLEL

    !!$ free temporaries
    deallocate(temps)
#undef D_UPW
#undef LEFT1
#undef RIGHT1
#undef LEFTPLUS
#undef LEFTMINUS
#undef RIGHTPLUS
#undef RIGHTMINUS
#undef CHECK
#undef MON
#undef STEEP
#undef VD_UPW
#undef VLEFT1
#undef VRIGHT1
#undef VLEFTPLUS
#undef VLEFTMINUS
#undef VRIGHTPLUS
#undef VRIGHTMINUS
#undef VCHECK
#undef VMON
#undef VSTEEP
#undef TSTOP
#undef TSTART

    return
  
  end subroutine SimplePPM_1d_dir

! tail of preprocessor loop, check only case in pre-processor loop

#if defined(SimplePPM_1d_z_dMHD_done)
end module PPM_SimplePPM_1d_XYZ
#endif

#undef IND3D
#undef VIND3D
#undef END_LOOPS
#undef BEGIN_LOOPS
#undef SimplePPM_1d_dir
#undef HAVE_BFIELD
#undef DIVCLEAN
#undef NTEMPS

! Include itself if we aren't at the end
#if defined(SimplePPM_1d_z_dMHD_done)
#  undef SimplePPM_1d_x_done
#  undef SimplePPM_1d_y_done
#  undef SimplePPM_1d_z_done
#  undef SimplePPM_1d_x_MHD_done
#  undef SimplePPM_1d_y_MHD_done
#  undef SimplePPM_1d_z_MHD_done
#  undef SimplePPM_1d_x_dMHD_done
#  undef SimplePPM_1d_y_dMHD_done
#  undef SimplePPM_1d_z_dMHD_done
#else
#include "/work2/09364/jakedoh/frontera/CactusMG/Cactus/arrangements/BHNSThorns/ReconstructPPM/src/PPM_SimplePPM_1d_XYZ.F90"
!! We are including ourselves again.
#endif
