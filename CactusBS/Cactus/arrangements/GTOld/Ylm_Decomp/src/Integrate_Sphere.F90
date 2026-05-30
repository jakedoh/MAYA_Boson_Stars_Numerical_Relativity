!  /*@@
!  @file      Integrate_Sphere.F90
!  @date      unknown
!  @author    unknown
!  @desc
!  @enddesc
!  @version   $Id: Integrate_Sphere.F90,v 1.30 2005/03/23 19:58:52 herrmann Exp $
!  @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


! function to precompute the spin-weight -2 Ylm array
! only called initially
subroutine Ylm_precompute_swm2(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: il,im,imarr,i,j, mmodes
  CCTK_REAL :: theta,phi,sint,cost

  CCTK_REAL,dimension(2) :: sp2Ylm,sm2Ylm,Ylm

  call CCTK_INFO("precomputing the spin-weight -2 array")
  swm2_ylm_pre=zero

  ! Loop over l-modes
  loop_l: do il = 2,l_mode

     ! Loop over m-modes
     loop_m: do im = -il,il
        imarr=il+im+1

        ! loop over points on surface
        loop_phi: do j = 1, maxnphi
           loop_theta: do i = 1, maxntheta
              
              theta=ctheta_max(i,j)
              phi=cphi_max(i,j)
              sint=sintheta_max(i,j)
              cost=costheta_max(i,j)
              
              call s2ylm_spher_combs(theta,phi,sint,cost,il,im,sp2Ylm,sm2Ylm,Ylm)
              
              swm2_ylm_pre(1,i,j,imarr,il)=sm2Ylm(1)
              swm2_ylm_pre(2,i,j,imarr,il)=sm2Ylm(2)

           end do loop_theta
        end do loop_phi
     end do loop_m
  end do loop_l

end subroutine Ylm_precompute_swm2





!/*@@
!  @routine    Ylm_Integrate_Sphere
!  @date       unknown
!  @author     unknown
!  @desc
!              Compute Regge Wheeler quantities and from them the Moncrief
!              Qeven, Qodd functions.
!  @enddesc
!@@*/
subroutine Ylm_Integrate_Sphere(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, il, im, ilarr, imarr
  CCTK_INT :: status, ierr
  CCTK_REAL :: st, ist, st2, ct2,ct

  CCTK_REAL,dimension(2) :: Ylm,Y1,Y2,Y3,Y4, &
                            h1,H2,K,G,c1,c2,dG,dK,dc2

  CCTK_REAL :: dtp

  CCTK_REAL,dimension(:),allocatable :: sum_var
  CCTK_REAL,dimension(:,:),allocatable :: output_tmp_2d
  CCTK_REAL,dimension(:,:,:),allocatable :: spherharm_2d
  CCTK_REAL :: sm2Yl2m0,Yl2m0, phi, theta, cost, sint
  CCTK_INT :: vind
  CCTK_REAL :: sth,cth,smn2Yl2m0,prefac
  CCTK_COMPLEX :: Yl2m2,Yl2mn2, smn2Yl2mn2,smn2Yl2m2, sm2Yl2mn2,sm2Yl2m2, tempc
  CCTK_COMPLEX :: cmplx_tmp
  CCTK_REAL,dimension(-2:2,0:2,-2:2) :: sylmro
  CCTK_COMPLEX,dimension(-2:2,0:2,-2:2) :: sylmco
  CCTK_REAL,dimension(-2:2) :: sylmr, sylmi
  CCTK_COMPLEX,dimension(-2:2) :: sylmc
  CCTK_REAL,dimension(2) :: sm2ylm, sp2ylm

  CCTK_REAL :: sqrval,xnorm,ynorm,drad, tmpval, dOmega

  CCTK_REAL :: rpart, ipart

  character(len=1024) :: infoline ! xxx len limited

  CCTK_REAL :: intweight, iwtheta, iwphi

  CCTK_INT :: skip_l, skip_m

! _________________________________________________________________


  if (verbose>4) &
    call CCTK_INFO("Compute Ylm integrals")

  if (do_nothing == 1) then
    return
  end if

  allocate(sum_var(nrdecompvars))
  sum_var=zero

  ! 2D tmp for output
  if (output_2d_surface>0) allocate(output_tmp_2d(cntheta,cnphi))
  if (output_2d_spherical_harmonics>0) allocate(spherharm_2d(cntheta,cnphi,6))

  dtp= dtheta*dphi

  ! l-mode setup
  if (CCTK_EQUALS(mode_type,"specific mode")) then
    l_min = l_mode ; l_max = l_mode
  else if (CCTK_EQUALS(mode_type,"all modes")) then
    if (lmin .eq. -1) then
       l_min=2
    else
       l_min=lmin
    end if
    if (lmax .eq. -1) then
       l_max=l_mode
    else
       l_max=lmax
    end if
  end if

  l_step=1
  m_step=1

  if (verbose > 3) then
    write(infoline,'(A29,I2,I2,I2)') 'mode setup: l [min,max,step]:',l_min,l_max,l_step
    call CCTK_INFO(infoline)
  end if

  if (.not.allocated(swm2_ylm_pre)) then
     call CCTK_WARN(1,"array swm2_ylm_pre not allocated!")
  end if

  ! Loop over l-modes
  loop_l: do il = l_min,l_max,l_step

    skip_l=0
    call Ylm_skip_l_mode_p(CCTK_PASS_FTOF,il,skip_l)
    if (skip_l>0) cycle

    call Ylm_m_mode_setup(CCTK_PASS_FTOF,il)

    ! Loop over m-modes
    loop_m: do im = m_min,m_max,m_step

      skip_m=0
      call Ylm_skip_m_mode_p(CCTK_PASS_FTOF,il,im,skip_m)
      if (skip_m>0) cycle

      sum_var=zero

      intweight=zero
      iwtheta=zero
      iwphi=zero

      ! loop over points on surface
      loop_phi1: do j = 1, cnphi
        ! integration weight
        if (j.eq.1 .or. j.eq.cnphi) then
           iwphi=3d0/8d0
        else if (j.eq.2 .or. j.eq.cnphi-1) then
           iwphi=7d0/6d0
        else if (j.eq.3 .or. j.eq.cnphi-2) then
           iwphi=23d0/24d0
        else
           iwphi=one
        end if

        loop_theta1: do i = 1, cntheta
          if (i.eq.1 .or. i.eq.cntheta) then
             iwtheta=zero
          else if (i.eq.2 .or. i.eq.cntheta-1) then
             iwtheta=55d0/24d0
          else if (i.eq.3 .or. i.eq.cntheta-2) then
             iwtheta=-1d0/6d0
          else if (i.eq.4 .or. i.eq.cntheta-3) then
             iwtheta=11d0/8d0
          else
             iwtheta=one
          end if

          intweight=iwphi*iwtheta

          theta=ctheta(i,j)
          phi=cphi(i,j)
          sint=sintheta(i,j)
          cost=costheta(i,j)

          ! populate ylms array, first index contains spin weight (either: -2,0,2)
          ! spin weight +/-2
          sylmr=zero
          sylmi=zero
          ! fuck fortran: cmplx() needs a kind qualifier otherwise
          ! it uses freaking single precision...
          sylmc=cmplx(zero,zero,8)
          if (use_precomputed_spwght==0) then
             call s2ylm_spher_combs(theta,phi,sint,cost,il,im,sp2Ylm,sm2Ylm,Ylm)
          else
             imarr=il+im+1
             sm2Ylm(1)=swm2_ylm_pre(1,i,j,imarr,il)
             sm2Ylm(2)=swm2_ylm_pre(2,i,j,imarr,il)
          end if
          ! FIXME : spin weight=1,-1 missing
          sylmr(-2)=sm2Ylm(1)
          sylmi(-2)=sm2Ylm(2)
          sylmr(0)=Ylm(1)
          sylmi(0)=Ylm(2)
          sylmr(2)=sp2Ylm(1)
          sylmi(2)=sp2Ylm(2)

          ! copy spherical harmonics for output if necessary
          if (output_2d_spherical_harmonics>0 .and. cctk_iteration .eq. 0) then
             spherharm_2d(i,j,1)=sm2Ylm(1)
             spherharm_2d(i,j,2)=sm2Ylm(2)
             spherharm_2d(i,j,3)=   Ylm(1)
             spherharm_2d(i,j,4)=   Ylm(2)
             spherharm_2d(i,j,5)=sp2Ylm(1)
             spherharm_2d(i,j,6)=sp2Ylm(2)
          end if

          ! fuck fortran: cmplx() needs a kind qualifier otherwise
          ! it uses freaking single precision...
          ! We want Ylm^* Psi_4 so here we take the cmplox conjugate
          sylmc=cmplx(sylmr,-sylmi,8)

          if (verbose>5) then
             print*,"compare the following"
             print*,"   for l,m=",il,im
             print*,"   for theta,phi=",theta,phi
             print*,"   sylmr(-2)...",sylmr(-2)
             print*,"   sylmr(0)...", sylmr(0)
             print*,"   sylmr(+2)...",sylmr(2)
             print*,"   sylmc(-2)...",sylmc(-2)
             print*,"   sylmc(0)...", sylmc(0)
             print*,"   sylmc(+2)...",sylmc(2)
          end if

          ! solid angle element times integration weight from simpson formula
          dOmega=sintheta(i,j)*dtp*intweight

          vind=1
          ! loop over variables
          do while (vind<=nrdecompvars)
            ! first check if we have a correct spin weight:
            if (      spwght(vind) .ne. 2 .and. spwght(vind) .ne. 0 &
                .and. spwght(vind) .ne. -2 ) then
               if (verbose>4) then
                  call CCTK_WARN(1,"wrong spin weight")
               end if
               cycle
            end if

            ! real gf
            if (cmplgf(vind) .eq. -1) then
              tmpval=interp_arrays(i,j, vind)*sylmr(spwght(vind))*dOmega
              tmpval=sylmr(spwght(vind))*dOmega
              sum_var(vind) = sum_var(vind)+tmpval
              integrand(i,j,vind)  = tmpval
              if (verbose>5) then
                 print*,"      [real] tmpv=",tmpval," sum_var(",vind,")=",sum_var(vind)
              end if

            ! complex gf
            ! imag part is in vind+1 entry
            else if (cmplgf(vind) .eq. 1) then
              if (vind==nrdecompvars) then
                print*,'FIXME : vind,nrdecompvars',vind,nrdecompvars
                call CCTK_WARN(1,"last variable flagged as complex. this can't be...")
              else
                ! fuck fortran: cmplx() needs a kind qualifier otherwise
                ! it uses freaking single precision...
                cmplx_tmp=cmplx(interp_arrays(i,j, vind), &
                     interp_arrays(i,j, vind+1),8)
                tempc=cmplx_tmp*sylmc(spwght(vind))
                tmpval=Real(tempc)*dOmega
                if (verbose>5) then
                   print*,"     [cmpl] real =",tmpval
                end if
                sum_var(vind) = sum_var(vind)+tmpval
                integrand(i,j,vind)  = tmpval
                tmpval=AImag(tempc)*dOmega
                if (verbose>5) then
                   print*,"     [cmpl] imag =",tmpval
                end if

                ! use next slot for imag part
                vind=vind+1
                sum_var(vind) = sum_var(vind)+tmpval
                integrand(i,j,vind)  = tmpval
                if (verbose>5) then
                   print*,"      [cmplx] sum_var(",vind-1,vind,")=",sum_var(vind-1),sum_var(vind)
                end if
              end if
            else
              call CCTK_WARN(1,"unknown var type, should be real or complex")
            end if
            vind=vind+1
          end do
        end do loop_theta1
      end do loop_phi1

      ! finally time for some IO
      ! now: sum_var() contains the value of the intelgrals
      ! output all detectors
      vind=1
      do while (vind<=nrdecompvars)
         if (verbose>3) then
            call CCTK_INFO("Output information for the following variable:")
            call CCTK_PrintVar(varindices(vind))
            write(infoline,'(A9,I2,A1,I2,A1)') '  (l,m)=(',il,',',im,')'
            call CCTK_INFO(infoline)
            write(infoline,*) '    integral = ', sum_var(vind)
            call CCTK_INFO(infoline)
         end if

         if (sum_var(vind).eq.zero.and.verbose>2) then
            call CCTK_INFO("Note: found exact zero in this variable:")
            call CCTK_PrintVar(varindices(vind))
         end if

         ! we might also be interested in the actual values of the
         ! 2D surface integrands over the sphere
         if (output_2d_surface>0) then
            ! copy the data into a temp
            output_tmp_2d=zero
            output_tmp_2d=integrand(:,:,vind)

            ! the integrand
            call write_surface_quantities_mode(ierr, cctk_time, &
                 ylm_files_exist, &
                 first_detector_this_time,&
                 current_detector,current_detector_radius, &
                 varindices(vind), il, im, &
                 output_tmp_2d(1,1), 1, cntheta, cnphi)
            if (ierr < 0) then
               call CCTK_WARN(1,"couldn't write ylm integrand")
            end if

            ! copy the data into a temp
            output_tmp_2d=zero
            output_tmp_2d=interp_arrays(:,:,vind)
            ! the interpolated function
            call write_surface_quantities_mode(ierr, cctk_time, &
                 ylm_files_exist, &
                 first_detector_this_time,&
                 current_detector,current_detector_radius, &
                 varindices(vind), il, im, &
                 output_tmp_2d(1,1), 0, cntheta, cnphi)
            if (ierr < 0) then
               call CCTK_WARN(1,"couldn't write ylm interpolated array")
            end if

         end if

         ! initially we can write out the spherical harmonics
         if (output_2d_spherical_harmonics>0 .and. cctk_iteration.eq.0) then
            ! spinw -2 real
            output_tmp_2d=zero 
            output_tmp_2d=spherharm_2d(:,:,1) 
            call write_surface_quantities_mode(ierr, cctk_time, & 
                 ylm_files_exist, & 
                 first_detector_this_time,& 
                 current_detector,current_detector_radius, & 
                 varindices(vind), il, im, & 
                 output_tmp_2d(1,1), 2, cntheta,cnphi) 
            if (ierr < 0) then 
               call CCTK_WARN(1,"couldn't write ylm spherical harmonics array") 
            end if

            ! spinw -2 imag
            output_tmp_2d=zero 
            output_tmp_2d=spherharm_2d(:,:,2) 
            call write_surface_quantities_mode(ierr, cctk_time, & 
                 ylm_files_exist, & 
                 first_detector_this_time,& 
                 current_detector,current_detector_radius, & 
                 varindices(vind), il, im, & 
                 output_tmp_2d(1,1), 3, cntheta, cnphi) 
            if (ierr < 0) then 
               call CCTK_WARN(1,"couldn't write ylm spherical harmonics array") 
            end if

            ! spinw 0 real
            output_tmp_2d=zero 
            output_tmp_2d=spherharm_2d(:,:,3) 
            call write_surface_quantities_mode(ierr, cctk_time, & 
                 ylm_files_exist, & 
                 first_detector_this_time,& 
                 current_detector,current_detector_radius, & 
                 varindices(vind), il, im, & 
                 output_tmp_2d(1,1), 4, cntheta,cnphi) 
            if (ierr < 0) then 
               call CCTK_WARN(1,"couldn't write ylm spherical harmonics array") 
            end if

            ! spinw 0 imag
            output_tmp_2d=zero 
            output_tmp_2d=spherharm_2d(:,:,4) 
            call write_surface_quantities_mode(ierr, cctk_time, & 
                 ylm_files_exist, & 
                 first_detector_this_time,& 
                 current_detector,current_detector_radius, & 
                 varindices(vind), il, im, & 
                 output_tmp_2d(1,1), 5,cntheta,cnphi) 
            if (ierr < 0) then 
               call CCTK_WARN(1,"couldn't write ylm spherical harmonics array") 
            end if

            ! spinw +2 real
            output_tmp_2d=zero 
            output_tmp_2d=spherharm_2d(:,:,5) 
            call write_surface_quantities_mode(ierr, cctk_time, & 
                 ylm_files_exist, & 
                 first_detector_this_time,& 
                 current_detector,current_detector_radius, & 
                 varindices(vind), il, im, & 
                 output_tmp_2d(1,1), 6,cntheta,cnphi) 
            if (ierr < 0) then 
               call CCTK_WARN(1,"couldn't write ylm spherical harmonics array") 
            end if

            ! spinw +2 imag
            output_tmp_2d=zero 
            output_tmp_2d=spherharm_2d(:,:,6) 
            call write_surface_quantities_mode(ierr, cctk_time, & 
                 ylm_files_exist, & 
                 first_detector_this_time,& 
                 current_detector,current_detector_radius, & 
                 varindices(vind), il, im, & 
                 output_tmp_2d(1,1), 7,cntheta,cnphi) 
            if (ierr < 0) then 
               call CCTK_WARN(1,"couldn't write ylm spherical harmonics array") 
            end if
         end if

         ! MAIN OUTPUT OF DATA:
         ! call output function to write data
         ! real gf
         if (cmplgf(vind) .eq. -1) then
            call write_ylm_mode(ierr, cctk_time, ylm_files_exist, &
                 0, &
                 first_detector_this_time,&
                 current_detector,current_detector_radius, &
                 varindices(vind), &
                 il, im,&
                 sum_var(vind),&
                 zero)
         ! cmplx gf
         else if (cmplgf(vind) .eq. 1) then
           if (vind==nrdecompvars) then
             print*,'FIXME : vind,nrdecompvars',vind,nrdecompvars
             call CCTK_WARN(1,"last variable flagged as complex. this can't be...")
           else
             rpart=sum_var(vind)
             ipart=sum_var(vind+1)

             call write_ylm_mode(ierr, cctk_time, ylm_files_exist, &
                  1, &
                  first_detector_this_time,&
                  current_detector,current_detector_radius, &
                  varindices(vind), &
                  il, im,&
                  sum_var(vind),sum_var(vind+1) )
             vind=vind+1
           end if
         end if

         if (ierr < 0) then
            call CCTK_WARN(1,"couldn't write ylm data")
         end if
         
         vind=vind+1
      end do
    end do loop_m
  end do loop_l

  if (output_2d_surface>0) deallocate(output_tmp_2d)
  if (output_2d_spherical_harmonics>0) deallocate(spherharm_2d)

  deallocate(sum_var)
end subroutine Ylm_Integrate_Sphere



! ------------------------------------------------------------------
!
! Calculate the (l,m) spherical harmonic at given angular
! coordinates. This number is in general complex, and
!
! Ylm = Ylm(1) + i Ylm(2)
!
! where
!
!           a    ( 2 l + 1 (l-|m|)! )                      i m phi
! Ylm = (-1) SQRT( ------- -------  ) P_l|m| (cos(theta)) e
!                (   4 Pi  (l+|m|)! )
!
! and where
!
! a = m/2 (sign(m)+1)
!
! ------------------------------------------------------------------
subroutine Ylm_spherical_harmonic(l,m,theta,phi,cost,Ylm)

  use Ylm_Constants

  implicit none

! Input variables
  CCTK_INT,intent(in) :: l,m
  CCTK_REAL,intent(in) :: theta,phi,cost

! Output variables
  CCTK_REAL,intent(out) :: Ylm(2)

! Local variables
  CCTK_INT :: i
  CCTK_REAL :: a,fac,Ylm_plgndr
! _________________________________________________________________

  fac = one
  do i = l-abs(m)+1,l+abs(m)
        fac = fac*dble(i)
  end do
  fac = one/fac

! a = (-one)**((m*ISIGN(m,1)/abs(m)+m)/2)*SQRT(dble(2*l+1)
! &    /four/Pi*fac)*Ylm_plgndr(l,abs(m),cos(theta))

  a = (-one)**max(m,0)*sqrt(dble(2*l+1)/ &
      four/Pi*fac)*Ylm_plgndr(l,abs(m),cost)
  Ylm(1) = a*cos(dble(m)*phi)
  Ylm(2) = a*sin(dble(m)*phi)
end subroutine Ylm_spherical_harmonic


! __________________________________________________________________
!
! FIXME : CHECK THE LICENSE OF THIS ROUTINE
! The numerical recipes license does only allow the owner of the book to use the
! routines on a single screen and single machine. I doubt that supercomputers
! count as single machine and certainly the user cannot use the code on her
! workstation as well as on a cluster. I would suggest replacing all of this
! nonsense by the gsl routine gsl_sf_legendre_sphPlm  
! http://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
! which can be used freely "in-house". Note that the gsl routine also includes
! the sqrt prefactors. Carpet requires the GSL anyway.
! From Numerical Recipes FIXME
!
!   Calculates the associated Legendre polynomial Plm(x).
!   Here m and l are integers satisfying 0 <= m <= l,
!   while x lies in the range -1 <= x <= 1
!
! __________________________________________________________________
function Ylm_plgndr(l,m,x)

  use Ylm_Constants

  implicit none

! Input variables
  CCTK_INT,INTENT(IN) :: l,m
  CCTK_REAL,INTENT(IN) :: x

! Output variables
  CCTK_REAL :: Ylm_plgndr

! Local Variables
  CCTK_INT :: i,ll
  CCTK_REAL :: pmm,somx2,fact,pmmp1,pll

! __________________________________________________________________

  pmm = one

  if (m.gt.0) then
    somx2=sqrt((one-x)*(one+x))
    fact=one
    do i=1,m
      pmm  = -pmm*fact*somx2
      fact = fact+two
    end do
  end if

  if (l.eq.m) then
    Ylm_plgndr = pmm
  else
    pmmp1 = x*(two*m+one)*pmm
    if (l.eq.m+1) then
      Ylm_plgndr=pmmp1
    else
      do ll=m+2,l
        pll = ( x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm )/dble(ll-m)
        pmm = pmmp1
        pmmp1 = pll
      end do
      Ylm_plgndr = pll
    end if
  end if
end function Ylm_plgndr


! __________________________________________________________________
!
! Calculates the various combinations of spherical harmonics needed
! for the extraction (all are complex):     
!
!   Y  = Ylm
!   Y1 = Ylm,theta
!   Y2 = Ylm,phi
!   Y3 = Ylm,theta,theta-cot theta Ylm,theta-Ylm,phi,phi/sin^2 theta
!   Y4 = Ylm,theta,phi-cot theta Ylm,phi
!
! The local variables Yplus is the spherical harmonic at (l+1,m)
! __________________________________________________________________
subroutine Ylm_spher_harm_combs(theta,phi,sint,cost,l,m,Y,Y1,Y2,Y3,Y4)

  use Ylm_Constants

  implicit none

! Input variables
  CCTK_INT,intent(in) :: l,m
  CCTK_REAL,intent(in) :: theta,phi,sint,cost

! Output variables
  CCTK_REAL,DIMENSION(2),intent(out) :: Y,Y1,Y2,Y3,Y4

! Local variables
  CCTK_INT :: i
  CCTK_REAL ::  Yplus(2),rl,rm,cot_theta
! __________________________________________________________________

  rl = dble(l)
  rm = dble(m)

  cot_theta = cost/sint

  call Ylm_spherical_harmonic(l+1,m,theta,phi,cost,Yplus)

  ! Find Y
  call Ylm_spherical_harmonic(l  ,m,theta,phi,cost,Y)

  ! Find Y1
  do i = 1,2
    Y1(i) = -(rl+one)*cot_theta*Y(i)+Yplus(i)/sint* &
            sqrt(((rl+one)**2-rm**2)*(rl+half)/(rl+one+half))
  end do

  ! Find Y2
  Y2(1) = -rm*Y(2)
  Y2(2) =  rm*Y(1)

  ! Find Y3
  do i = 1,2
    Y3(i) = -two*cot_theta*Y1(i)+(two*rm*rm/(sint**2) &
            -rl*(rl+one))*Y(i)
  end do

  ! Find Y4
  Y4(1) = rm*(cot_theta*Y(2)-Y1(2))
  Y4(2) = rm*(Y1(1)-cot_theta*Y(1))
end subroutine Ylm_spher_harm_combs


subroutine s2ylm_spher_combs(theta,phi,sint,cost,l,m,sp2Ylm,sm2Ylm,Ylm)
  use Ylm_Constants

  implicit none

! Input variables
  CCTK_INT,intent(in) :: l,m
  CCTK_REAL,intent(in) :: theta,phi,sint,cost

! Output variables
  CCTK_REAL,DIMENSION(2),intent(out) :: sp2Ylm,sm2Ylm,Ylm

! Local variables
  CCTK_INT :: i
  CCTK_REAL,DIMENSION(2) :: Y,Y1,Y2,Y3,Y4
  CCTK_COMPLEX :: Wlm,Xlm
  CCTK_REAL :: fac

  call Ylm_spher_harm_combs(theta,phi,sint,cost,l,m,Y,Y1,Y2,Y3,Y4)

  ! fuck fortran: cmplx() needs a kind qualifier otherwise
  ! it uses freaking single precision...
  Wlm=cmplx(Y3(1),Y3(2),8)
  Xlm=two*cmplx(Y4(1),Y4(2),8)

  fac = one
  do i = l-1,l+2
        fac = fac*dble(i)
  end do
  fac = sqrt(one/fac)

  sp2Ylm(1)=fac*Real(Wlm+iunit/sint*Xlm)
  sp2Ylm(2)=fac*AImag(Wlm+iunit/sint*Xlm)

  sm2Ylm(1)=fac*Real(Wlm-iunit/sint*Xlm)
  sm2Ylm(2)=fac*AImag(Wlm-iunit/sint*Xlm)

  Ylm=Y

end subroutine s2ylm_spher_combs

subroutine Ylm_skip_l_mode_p(CCTK_ARGUMENTS,l,skip_l)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, intent(in) :: l
  CCTK_INT, intent(out) :: skip_l

  skip_l=0
  if (CCTK_EQUALS(l_even_odd,"all")) then
    skip_l=0
  else if (CCTK_EQUALS(l_even_odd,"even")) then
    skip_l=0
    if (mod(l+1,2).eq.0) skip_l=1
  else if (CCTK_EQUALS(l_even_odd,"odd")) then
    skip_l=0
    if (mod(l,2).eq.0) skip_l=1
  end if
end subroutine Ylm_skip_l_mode_p

subroutine Ylm_skip_m_mode_p(CCTK_ARGUMENTS,il,m,skip_m)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, intent(out) :: skip_m
  CCTK_INT,intent(in) :: m,il

  skip_m=0
  
  if (m < -il) then
     skip_m=1
     return
  end if
  if (m > il) then
     skip_m=1
     return
  end if

  if (CCTK_EQUALS(m_even_odd,"all")) then
    skip_m=0
  else if (CCTK_EQUALS(m_even_odd,"even")) then
    skip_m=0
    if (mod(m+1,2).eq.0) skip_m=1
  else if (CCTK_EQUALS(m_even_odd,"odd")) then
    skip_m=0
    if (mod(m,2).eq.0) skip_m=0
  end if
end subroutine Ylm_skip_m_mode_p

subroutine Ylm_m_mode_setup(CCTK_ARGUMENTS,il)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, intent(in) :: il

  character(len=1024) :: infoline ! xxx len limited

  ! m-mode setup (depends on l_mode, ie m_mode <= l_mode)
  if (CCTK_EQUALS(mode_type,"specific mode")) then
     m_min = m_mode ; m_max = m_mode
  else if(CCTK_EQUALS(mode_type,"all modes")) then
     if (mmin .eq. -2048) then ! ugly default value of -2048 xxx
        m_min = -m_mode 
     else
        m_min = mmin
     end if
     if (mmax .eq. -2048) then
        m_max = m_mode
     else
        m_max = mmax
     end if
     if (m_max>il) then
        m_max=il
     end if
     if (m_min< -il) then
        m_min=-il
     end if
  end if

  if (verbose > 3) then
     write(infoline,'(A29,I2,I2,I2)') '            m [min,max,step]:',m_min,m_max,m_step
     call CCTK_INFO(infoline)
  end if

end subroutine Ylm_m_mode_setup
