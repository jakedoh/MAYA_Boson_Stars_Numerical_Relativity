!/*@@
!  @file      MoncriefQ.F90
!  @date      unknown
!  @author    unknown
!  @desc
!             Compute Regge Wheeler quantities and from them the Moncrief
!             Qeven, Qodd functions.
!  @enddesc
!  @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

!/*@@
!  @routine    WavExtrL_MoncriefQ
!  @date       unknown
!  @author     unknown
!  @desc
!              Compute Regge Wheeler quantities and from them the Moncrief
!              Qeven, Qodd functions.
!  @enddesc
!@@*/
subroutine WavExtrL_MoncriefQ(CCTK_ARGUMENTS)

  use WavExtrLConstants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, il, im
  CCTK_INT :: status, istat
  integer ::  ierr, sumhandle
  CCTK_REAL :: st, ist
  CCTK_INT,dimension(2) :: lsh

  CCTK_REAL :: fac_h1, fac_H2, fac_G, fac_K, &
               fac_c1, fac_c2, fac_dG, fac_dK, fac_dc2

  CCTK_REAL :: gttcomb, gtpcomb, gppcomb

  CCTK_REAL :: Lambda

  CCTK_REAL,dimension(2) :: Ylm,Y1,Y2,Y3,Y4, &
                            h1,H2,K,G,c1,c2,dG,dK,dc2, &
                            ih1,iH2,iK,iG,ic1,ic2,idG,idK,idc2

  CCTK_REAL :: dtheta, dphi, dtp

  integer :: num_out_vals, num_in_fields, minus_one
  CCTK_REAL,dimension(18) :: out_vals, local_reduced_vals(18)

  CCTK_REAL :: lapse,cor_angle
  CCTK_INT :: marr


  character(len=80) :: infoline

! _________________________________________________________________

  if (verbose>4) &
    call CCTK_INFO("Compute Regge Wheeler quantities and from them Qeven, Qodd")

  if (do_nothing == 1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,my_out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
  end if

  if (calc_when_necessary .eq. 1) then
    if (cctk_time .lt. current_detector_radius-50) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    endif
    call CCTK_IsFunctionAliased(istat, "MergerHandler_WeHaveMerger")
    if (istat .eq. 1) then
      if (MergerHandler_WeHaveMerger() .eq. 1) then
        if (cctk_time .gt. MergerHandler_MergerTime()+current_detector_radius+ringdown_margin) then
          if (verbose>2) call CCTK_INFO("No time for this detector")
          return
        endif
      endif
    endif
  end if

  ! local shape of grid arrays on sphere
  call CCTK_GrouplshGN(status, cctkGH, 2, lsh, "WaveExtractL::surface_arrays")
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, "cannot get local size for surface arrays" )
  end if

  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon .ne. 0) then
    dphi = two*pi
  end if

  dtp= dtheta*dphi

  ! l-mode setup
  if (CCTK_EQUALS(mode_type,"specific mode")) then
    l_min = l_mode ; l_max = l_mode
  else if (CCTK_EQUALS(mode_type,"all modes")) then
    l_min = 2 ; l_max = l_mode
  end if


  if (cartoon .ne. 0) then
    l_step =2
    m_step =2
  else if (CCTK_EQUALS(domain,"full")) then
    l_step =1
    m_step =1
  else if (CCTK_EQUALS(domain,"octant")) then
    l_step =2
    m_step =2
  ! FIXME : look at symmetries
  else if (CCTK_EQUALS(domain,"quadrant")) then
    l_step =1
    m_step =1
  ! FIXME : REALLY ALL MODES?? Does the bitant symmetry kill some, ie odd modes ???
  else if (CCTK_EQUALS(domain,"bitant")) then
    l_step =1
    m_step =1
  end if

  if (verbose > 3) then
    write(infoline,'(A29,I2,I2,I2)') 'mode setup: l [min,max,step]:',l_min,l_max,l_step
    call CCTK_INFO(infoline)
  end if

  ! Loop over l-modes
  loop_l: do il = l_min,l_max,l_step

    ! m-mode setup (depends on l_mode, ie m_mode <= l_mode)
    if (CCTK_EQUALS(mode_type,"specific mode")) then
      m_min = -m_mode ; m_max = m_mode
    else if(CCTK_EQUALS(mode_type,"all modes")) then
      m_min = 0 ; m_max = m_mode
      if (m_mode>il) then
        m_max=il
      end if
    end if

    ! Factors independent of angular coordinates and of m-mode, but not l_mode
    fac_h1  = dri_drsch/dble(il*(il+1))
    fac_H2  = S_factor*dri_drsch**2
    fac_G   = one/(rsch**2*dble(il*(il+1)*(il-1)*(il+2)))
    fac_K   = half/rsch**2
    fac_c1  = dri_drsch/dble(il*(il+1))
    fac_c2  = two/dble(il*(il+1)*(il-1)*(il+2))
    fac_dG  = fac_G
    fac_dK  = fac_K
    fac_dc2 = fac_c2

    call CCTK_TimerStart(ierr,"MoncriefQ")

    if (verbose > 3) then
      write(infoline,'(A29,I2,I2,I2)') '            m [min,max,step]:',m_min,m_max,m_step
      call CCTK_INFO(infoline)
    end if

    ! Get sum reduction operator handle
    call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
    if ( sum_handle .lt. 0 ) then
      call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
    end if

    ! Loop over m-modes
    loop_m: do im = m_min,m_max,m_step
      ! Real parts of the Regge-Wheeler variables
      loop_phi1: do j = 1, lsh(2)
        loop_theta1: do i = 1, lsh(1)

          st  = sintheta(i,j)
          ist = one/st

          call WavExtrL_spher_harm_combs(ctheta(i,j),cphi(i,j),il,im,Ylm,Y1,Y2,Y3,Y4)

          h1i(i,j) =  st*grt(i,j)*Y1(1)+ist*grp(i,j)*Y2(1)
          H2i(i,j) =  st*grr(i,j)*Ylm(1)
          Gi(i,j)  =  (st*gtt(i,j)-ist*gpp(i,j))*Y3(1) &
                         +four*ist*gtp(i,j)*Y4(1)
          Ki(i,j)  =  (st*gtt(i,j)+ist*gpp(i,j))*Ylm(1)
          c1i(i,j) =  grp(i,j)*Y1(1)-grt(i,j)*Y2(1)
          c2i(i,j) =  (gtt(i,j)-ist**2*gpp(i,j))*Y4(1) &
                         -gtp(i,j)*Y3(1)

          gttcomb = dri_drsch*dr_gtt(i,j)-two/rsch*gtt(i,j)
          gtpcomb = dri_drsch*dr_gtp(i,j)-two/rsch*gtp(i,j)
          gppcomb = dri_drsch*dr_gpp(i,j)-two/rsch*gpp(i,j)

          dGi(i,j)  = (st*gttcomb-ist*gppcomb)*Y3(1) &
                        +four*ist*gtpcomb*Y4(1)
          dKi(i,j)  =  (st*gttcomb+ist*gppcomb)*Ylm(1)
          dc2i(i,j) =  (dr_gtt(i,j)-ist**2*dr_gpp(i,j))*Y4(1) &
                          -dr_gtp(i,j)*Y3(1)

          ! for m!=0 case we have imaginary part as well from the Ylm's
          if (im.ne.0) then
            ! switch signs - stupid convention
            Ylm(2)=-Ylm(2)
            Y1(2)=-Y1(2)
            Y2(2)=-Y2(2)
            Y3(2)=-Y3(2)
            Y4(2)=-Y4(2)
            ih1i(i,j) =  st*grt(i,j)*Y1(2)+ist*grp(i,j)*Y2(2)
            iH2i(i,j) =  st*grr(i,j)*Ylm(2)
            iGi(i,j)  =  (st*gtt(i,j)-ist*gpp(i,j))*Y3(2) &
                           +four*ist*gtp(i,j)*Y4(2)
            iKi(i,j)  =  (st*gtt(i,j)+ist*gpp(i,j))*Ylm(2)
            ic1i(i,j) =  grp(i,j)*Y1(2)-grt(i,j)*Y2(2)
            ic2i(i,j) =  (gtt(i,j)-ist**2*gpp(i,j))*Y4(2) &
                           -gtp(i,j)*Y3(2)

            idGi(i,j)  = (st*gttcomb-ist*gppcomb)*Y3(2) &
                          +four*ist*gtpcomb*Y4(2)
            idKi(i,j)  =  (st*gttcomb+ist*gppcomb)*Ylm(2)
            idc2i(i,j) =  (dr_gtt(i,j)-ist**2*dr_gpp(i,j))*Y4(2) &
                            -dr_gtp(i,j)*Y3(2)
          end if
        end do loop_theta1
      end do loop_phi1

      ! Integrations over the 2-sphere
      ! Note the abscence of sintheta which is already included
      ! in the above expressions
      int_tmp1=sym_factor*weights*dtp* h1i
      int_tmp2=sym_factor*weights*dtp* H2i
      int_tmp3=sym_factor*weights*dtp* Gi
      int_tmp4=sym_factor*weights*dtp* Ki
      int_tmp5=sym_factor*weights*dtp* c1i
      int_tmp6=sym_factor*weights*dtp* c2i
      int_tmp7=sym_factor*weights*dtp* dGi
      int_tmp8=sym_factor*weights*dtp* dKi
      int_tmp9=sym_factor*weights*dtp* dc2i

      local_reduced_vals(1) = sum(int_tmp1,weights.gt.1.0e-15)
      local_reduced_vals(2) = sum(int_tmp2,weights.gt.1.0e-15)
      local_reduced_vals(3) = sum(int_tmp3,weights.gt.1.0e-15)
      local_reduced_vals(4) = sum(int_tmp4,weights.gt.1.0e-15)
      local_reduced_vals(5) = sum(int_tmp5,weights.gt.1.0e-15)
      local_reduced_vals(6) = sum(int_tmp6,weights.gt.1.0e-15)
      local_reduced_vals(7) = sum(int_tmp7,weights.gt.1.0e-15)
      local_reduced_vals(8) = sum(int_tmp8,weights.gt.1.0e-15)
      local_reduced_vals(9) = sum(int_tmp9,weights.gt.1.0e-15)

      num_out_vals=1
      minus_one = -1
      sumhandle = sum_handle ! i.e., convert from CCTK_INT to integer

      if (im.eq.0) then
        num_in_fields=9

      else
        num_in_fields=18

        int_tmp10=sym_factor*weights*dtp* ih1i
        int_tmp11=sym_factor*weights*dtp* iH2i
        int_tmp12=sym_factor*weights*dtp* iGi
        int_tmp13=sym_factor*weights*dtp* iKi
        int_tmp14=sym_factor*weights*dtp* ic1i
        int_tmp15=sym_factor*weights*dtp* ic2i
        int_tmp16=sym_factor*weights*dtp* idGi
        int_tmp17=sym_factor*weights*dtp* idKi
        int_tmp18=sym_factor*weights*dtp* idc2i

        local_reduced_vals(10) = sum(int_tmp10,weights.gt.1.0e-15)
        local_reduced_vals(11) = sum(int_tmp11,weights.gt.1.0e-15)
        local_reduced_vals(12) = sum(int_tmp12,weights.gt.1.0e-15)
        local_reduced_vals(13) = sum(int_tmp13,weights.gt.1.0e-15)
        local_reduced_vals(14) = sum(int_tmp14,weights.gt.1.0e-15)
        local_reduced_vals(15) = sum(int_tmp15,weights.gt.1.0e-15)
        local_reduced_vals(16) = sum(int_tmp16,weights.gt.1.0e-15)
        local_reduced_vals(17) = sum(int_tmp17,weights.gt.1.0e-15)
        local_reduced_vals(18) = sum(int_tmp18,weights.gt.1.0e-15)

      end if

      call CCTK_ReduceLocArrayToArray1D(ierr, cctkGH, minus_one,&
                        sumhandle, local_reduced_vals(1:num_in_fields),&
                        out_vals(1:num_in_fields), num_in_fields,&
                        CCTK_VARIABLE_REAL)

      if (ierr.ne.0) then
        call CCTK_WARN(1,"The reduction of the MoncriefQ integrands failed!")
      end if

      h1(1)  = out_vals(1)
      H2(1)  = out_vals(2)
      G(1)   = out_vals(3)
      K(1)   = out_vals(4)
      c1(1)  = out_vals(5)
      c2(1)  = out_vals(6)
      dG(1)  = out_vals(7)
      dK(1)  = out_vals(8)
      dc2(1) = out_vals(9)

      h1(1)  = fac_h1 * h1(1)
      H2(1)  = fac_h2 * H2(1)
      G(1)   = fac_G  * G(1)
      K(1)   = fac_K  * K(1)  +dble(il*(il+1))*half*G(1)
      c1(1)  = fac_c1 * c1(1)
      c2(1)  = fac_c2 * c2(1)
      dG(1)  = fac_dG * dG(1)
      dK(1)  = fac_dK * dK(1) +dble(il*(il+1))*half*dG(1)
      dc2(1) = fac_dc2* dc2(1)


      if (verbose >4) then
        print*,'Real Quantities'
        print*,'h1',h1(1)
        print*,'H2',H2(1)
        print*,'G',G(1)
        print*,'K',K(1)
        print*,'c1',c1(1)
        print*,'c2',c2(1)
        print*,'dG',dG(1)
        print*,'dK',dK(1)
        print*,'dc2',dc2(1)
      end if



      if (im.ne.0) then
        h1(2)  = out_vals(10)
        H2(2)  = out_vals(11)
        G(2)   = out_vals(12)
        K(2)   = out_vals(13)
        c1(2)  = out_vals(14)
        c2(2)  = out_vals(15)
        dG(2)  = out_vals(16)
        dK(2)  = out_vals(17)
        dc2(2) = out_vals(18)

        h1(2)  = fac_h1 * h1(2)
        H2(2)  = fac_h2 * H2(2)
        G(2)   = fac_G  * G(2)
        K(2)   = fac_K  * K(2)  +dble(il*(il+1))*half*G(2)
        c1(2)  = fac_c1 * c1(2)
        c2(2)  = fac_c2 * c2(2)
        dG(2)  = fac_dG * dG(2)
        dK(2)  = fac_dK * dK(2) +dble(il*(il+1))*half*dG(2)
        dc2(2) = fac_dc2* dc2(2)
      end if

      ! Moncrief Q
      Lambda = dble((il-1)*(il+2))+three*(one-S_factor)

      ! m index into array : fortran: 1,2,3 index -> -m_max,-m_max+1,...
      marr=-m_min+im +1
      if (il<1 .or. il>l_mode .or. marr<1 .or. marr>m_mode+1) call CCTK_WARN (0, "internal error")
      Qodd_Re(il,marr) = sqrt(two*dble((il+2)*(il+1)*il*(il-1)))* &
                      S_factor/rsch*(c1(1)+half*(dc2(1)-two/rsch*c2(1)))


      Qeven_Re(il,marr) = one/Lambda*sqrt(two*dble((il-1)*(il+2))/ &
                          dble(il*(il+1)))*( dble(il*(il+1))*S_factor* &
              (rsch**2*dG(1)-two*h1(1))+two*rsch*S_factor*(H2(1)-rsch*dK(1)) &
                             +Lambda*rsch*K(1) )

      if (im.ne.0) then
        Qodd_Im(il,marr) = sqrt(two*dble((il+2)*(il+1)*il*(il-1)))* &
                        S_factor/rsch*(c1(2)+half*(dc2(2)-two/rsch*c2(2)))


        Qeven_Im(il,marr) = one/Lambda*sqrt(two*dble((il-1)*(il+2))/ &
                            dble(il*(il+1)))*( dble(il*(il+1))*S_factor* &
                (rsch**2*dG(2)-two*h1(2))+two*rsch*S_factor*(H2(2)-rsch*dK(2)) &
                               +Lambda*rsch*K(2) )
      end if

      if (verbose>3) then
        write(infoline,'(A9,I2,A1,I2,A1)') '  (l,m)=(',il,',',im,')'
        call CCTK_INFO(infoline)
        write(infoline,'(A10,G20.8,A7,G20.8)') '    Qeven=', &
                         Qeven_Re(il,marr),', Qodd=',Qodd_Re(il,marr)
        call CCTK_INFO(infoline)
        if (im.ne.0) then
          write(infoline,'(A16,G20.8,A7,G20.8)') '    imag: Qeven=', &
                         Qeven_Im(il,marr),', Qodd=',Qodd_Im(il,marr)
          call CCTK_INFO(infoline)
        end if
      end if
      if (current_detector<1 .or. current_detector>maximum_detector_number .or. il<1 .or. il>l_mode .or. marr<1 .or. marr>2*m_mode+1) call CCTK_WARN (0, "internal error")
      Qodd_Re_Array(current_detector,il,marr) = Qodd_Re(il,marr)
      Qeven_Re_Array(current_detector,il,marr) = Qeven_Re(il,marr)

      if (im.ne.0) then
        Qodd_Im_Array(current_detector,il,marr) = Qodd_Im(il,marr)
        Qeven_Im_Array(current_detector,il,marr) = Qeven_Im(il,marr)
      else
        Qodd_Im_Array(current_detector,il,marr) = zero
        Qeven_Im_Array(current_detector,il,marr) = zero
      end if
    end do loop_m
  end do loop_l

  call CCTK_TimerStop(ierr,"MoncriefQ")

  if (phicorotate .ne. 0) then
    call CCTK_INFO("reset corotation on phi at the end of MoncriefQ for next step")
    lapse=one
    cor_angle=cctk_time*rotation_omega/lapse

!   undo rotation transformation on phi itself.
    do j = 1, lsh(2)
      cphi(:,j) = cphi(:,j) + cor_angle
    end do

    sinphi = sin(cphi)
    cosphi = cos(cphi)
  end if

end subroutine WavExtrL_MoncriefQ



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
subroutine WavExtrL_spherical_harmonic(l,m,theta,phi,Ylm)

  use WavExtrLConstants

  implicit none

! Input variables
  CCTK_INT :: l,m
  CCTK_REAL :: theta,phi

! Output variables
  CCTK_REAL :: Ylm(2)

! Local variables
  CCTK_INT :: i
  CCTK_REAL :: a,fac,WavExtrL_plgndr
! _________________________________________________________________

  fac = one
  do i = l-abs(m)+1,l+abs(m)
        fac = fac*dble(i)
  end do
  fac = one/fac

! a = (-one)**((m*ISIGN(m,1)/abs(m)+m)/2)*SQRT(dble(2*l+1)
! &    /four/Pi*fac)*WavExtrL_plgndr(l,abs(m),cos(theta))

  a = (-one)**max(m,0)*sqrt(dble(2*l+1)/ &
      four/Pi*fac)*WavExtrL_plgndr(l,abs(m),cos(theta))
  Ylm(1) = a*cos(dble(m)*phi)
  Ylm(2) = a*sin(dble(m)*phi)
end subroutine WavExtrL_spherical_harmonic


! __________________________________________________________________
!
! FIXME : CHECK THE LICENSE OF THIS ROUTINE
! From Numerical Recipes FIXME
!
!   Calculates the associated Legendre polynomial Plm(x).
!   Here m and l are integers satisfying 0 <= m <= l,
!   while x lies in the range -1 <= x <= 1
!
! __________________________________________________________________
function WavExtrL_plgndr(l,m,x)

  use WavExtrLConstants

  implicit none

! Input variables
  CCTK_INT,INTENT(IN) :: l,m
  CCTK_REAL,INTENT(IN) :: x

! Output variables
  CCTK_REAL :: WavExtrL_plgndr

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
    WavExtrL_plgndr = pmm
  else
    pmmp1 = x*(two*m+one)*pmm
    if (l.eq.m+1) then
      WavExtrL_plgndr=pmmp1
    else
      do ll=m+2,l
        pll = ( x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm )/dble(ll-m)
        pmm = pmmp1
        pmmp1 = pll
      end do
      WavExtrL_plgndr = pll
    end if
  end if
end function WavExtrL_plgndr


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
!
! All the return values are 2 dim to account for real and complex part.
! __________________________________________________________________
subroutine WavExtrL_spher_harm_combs(theta,phi,l,m,Y,Y1,Y2,Y3,Y4)

  use WavExtrLConstants

  implicit none

! Input variables
  CCTK_INT :: l,m
  CCTK_REAL :: theta,phi

! Output variables
  CCTK_REAL,DIMENSION(2) :: Y,Y1,Y2,Y3,Y4

! Local variables
  CCTK_INT :: i
  CCTK_REAL ::  Yplus(2),rl,rm,cot_theta
! __________________________________________________________________

  rl = dble(l)
  rm = dble(m)

  cot_theta = cos(theta)/sin(theta)

  call WavExtrL_spherical_harmonic(l+1,m,theta,phi,Yplus)

  ! Find Y
  call WavExtrL_spherical_harmonic(l,m,theta,phi,Y)

  ! Find Y1
  do i = 1,2
    Y1(i) = -(rl+one)*cot_theta*Y(i)+Yplus(i)/sin(theta)* &
            sqrt(((rl+one)**2-rm**2)*(rl+half)/(rl+one+half))
  end do

  ! Find Y2
  Y2(1) = -rm*Y(2)
  Y2(2) =  rm*Y(1)

  ! Find Y3
  do i = 1,2
    Y3(i) = -two*cot_theta*Y1(i)+(two*rm*rm/(sin(theta)**2) &
            -rl*(rl+one))*Y(i)
  end do

  ! Find Y4
  Y4(1) = rm*(cot_theta*Y(2)-Y1(2))
  Y4(2) = rm*(Y1(1)-cot_theta*Y(1))
end subroutine WavExtrL_spher_harm_combs
