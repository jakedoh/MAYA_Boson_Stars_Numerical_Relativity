/*@@
c  @file      RotatingDBH.F
c  @date      8 June 1999
c  @author    Ryoji Takahashi
c  @desc 
c  
c  @enddesc 
c@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

c/*@@
c  @routine    RotatingDBHIVP
c  @date       
c  @author     
c  @desc 
c  
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c@@*/

      subroutine RotatingDBHIVP(CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      real*8 deta,dq,dphi
      real*8, allocatable :: ac(:,:,:),ae(:,:,:),aw(:,:,:),an(:,:,:),
     $     as(:,:,:),aq(:,:,:),ab(:,:,:),rhs(:,:,:),Ksq(:,:,:),
     $     psisph(:,:,:),psip(:,:,:),psipp(:,:,:),detapsisph(:,:,:),
     $     dqpsisph(:,:,:),dphipsisph(:,:,:),detaetapsisph(:,:,:),
     $     detaqpsisph(:,:,:),detaphipsisph(:,:,:),dqqpsisph(:,:,:),
     $     dqphipsisph(:,:,:),dphiphipsisph(:,:,:)
      real*8, allocatable :: etagrd(:),qgrd(:),phigrd(:)
      real*8 o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,o17,
     $     o18,o19,o20,o21,o22,o23,o24,o25,o26,o27,o28,o29,o30,o31,o32,
     $     o33,o34,o35,o36,o37,o38,o39,o40,o41,o42,o43,o44,o45,o46,o47,
     $     o48,o49,o50,o51,o52,o53,o54,o55,o56,o57,o58,o59,o60,o61,o62,
     $     o63,o64,o65,o66,o67,o68,o69,o70,o71,o72,o73,o74,o75,o76,o77,
     $     o78,o79,o80,o81,o82,o83,o84,o85,o86,o87,o88,o89,o90,o91,o92,
     $     o93,o94,o95,o96,o97,o98,o99,o100,o101,o102,o103,o104,o105,
     $     o106,o107,o108,o109,o110,o111,o112,o113,o114,o115,o116,o117,
     $     o118,o119,o120,o121,o122,o123,o124
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,
     $     t11,t12,t13,t14,t15,t16
      real*8 gtil,dngtil,dnngtil,dnnngtil,dnnnngtil,dnnnnngtil,
     $     dnnnnnngtil,dnnnnnnngtil
      real*8 rhsmax,rmax,adm
      real*8,parameter :: rbh_tol = 1.0d-7,rbh_eps = 1.0d-10
      integer,parameter :: itmax = 30 
      real*8 pi
      integer :: ne,nq,np
      integer :: nx,ny,nz
      integer i,j,k,it,ier,nquads,nocts,ntries
      integer npoints,ierror
      integer make_conformal_derivs

      integer       param_table_handle, interp_handle
      character(30) options_string
      CCTK_REAL,    dimension(3)  :: coord_origin, coord_delta
      CCTK_POINTER, dimension(3)  :: interp_coords
      CCTK_INT,     dimension(3 ) :: in_array_dims
      CCTK_POINTER, dimension(10) :: in_arrays, out_arrays
      CCTK_INT,     dimension(10), parameter :: type_codes = CCTK_VARIABLE_REAL


c Check if we should create and store conformal factor stuff 

      if (CCTK_EQUALS(metric_type, "static conformal")) then
    
         conformal_state = 1

         if(CCTK_EQUALS(conformal_storage,"factor+derivs")) then

            conformal_state = 2
            make_conformal_derivs = 1

          else if (CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs")) then

            conformal_state = 3
            make_conformal_derivs = 1

          endif

       else

         make_conformal_derivs = 0

       endif

      pi = 4.0d0*atan(1.0d0)

c     DO NOT use integer*4
      ne = neta
      nq = ntheta
      np = nphi

c     Set up the grid spacings
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)      
       
c     Distorted Rotating BH parameters (for Extrinsic curveture!)
c     
      print *,"Brill wave + Rotating Distorted BH solve"
      write(*,123)amp,eta0,sigma,byJ,mm
      print*,'etamax=',etamax
 123  format(1x, 'Pars: amp',f8.5,' eta0',f8.5,' sigma',f8.5,' byJ',f9.5
     $     ,' mm',f8.5)
      
c     Sovle on this sized cartesian grid
c     3D grid size NE x NT
c     Add 2 zones for eta coordinate and 4 for theta
c     and phi coordenate.
      ne = ne+2
      nq = nq+2
      np = np+2
c     
      allocate(ac(ne,nq,np),ae(ne,nq,np),aw(ne,nq,np),an(ne,nq,np),
     $     as(ne,nq,np),aq(ne,nq,np),ab(ne,nq,np),rhs(ne,nq,np),
     $     Ksq(ne,nq,np),psisph(ne,nq,np),psip(ne,nq,np),
     $     psipp(ne,nq,np),detapsisph(ne,nq,np),dqpsisph(ne,nq,np),
     $     dphipsisph(ne,nq,np),detaetapsisph(ne,nq,np),
     $     detaqpsisph(ne,nq,np),detaphipsisph(ne,nq,np),
     $     dqqpsisph(ne,nq,np),dqphipsisph(ne,nq,np),
     $     dphiphipsisph(ne,nq,np))
      allocate(etagrd(ne),qgrd(nq),phigrd(np))
c     
c     Initialize some array
c     
      nocts = 4
      nquads = 2      
      dphi = nocts*0.5*pi/(np-2)
      dq = nquads*0.5*pi/(nq-2)
      deta = etamax/(ne-3)
      
      do k = 1,np
         phigrd(k) = (k-1.5)*dphi
      enddo
      do j = 1,nq
         qgrd(j) = (j-1.5)*dq
      enddo
      do i=1,ne
         etagrd(i) = (i-2)*deta
      enddo
      
c     
c     Specify the initial guess for the conformal factor:
c     
      psip = 0.
      psipp = 0.
      do k = 1,np
         do j = 1,nq
            do i = 1,ne
               psisph(i,j,k) = 2.*cosh(0.5*etagrd(i))
            enddo
         enddo
      enddo
c     
c     Compute K^2:
c     
      do k = 1,np
         do j = 1,nq
            do i = 1,ne
#include "gauss.x"
#include "ksq_odd.x"
            enddo
         enddo
      enddo
c     
c     Solve the Hamiltonian constraint for the conformal factor:
c     
       ntries = 0
       
      do it = 1,itmax+1
         
         ac = 0.
         ae = 0.
         aw = 0.
         an = 0.
         as = 0.
         aq = 0.
         ab = 0.
         rhs = 0.
                  
c     
         do k = 2,np-1
            do j = 2,nq-1
               do i = 2,ne-1
                  rhs(i,j,k) =
     $                 (psip(i+1,j,k)-2.*psip(i,j,k)+psip(i-1,j,k))/
     $                 deta**2+(psip(i,j+1,k)-2.*psip(i,j,k)+
     $                 psip(i,j-1,k))/dq**2+0.5*(psip(i,j+1,k)-
     $                 psip(i,j-1,k))/(dq*tan(qgrd(j)))+
     $                 (psip(i,j,k+1)-2.*psip(i,j,k)+psip(i,j,k-1))/
     $                 (dphi**2*sin(qgrd(j))**2)-0.25*psip(i,j,k)+0.125*
     $                 Ksq(i,j,k)/(psisph(i,j,k)+psip(i,j,k))**7
               enddo
            enddo
         enddo
c     
         rhs = -rhs
         
c     
c     See if the residual is small enough:
c     
         rhsmax = 0.
c     
         do k = 2,np-1
            do j = 2,nq-1
               do i = 2,ne-1
                  rhsmax = max(rhsmax,abs(rhs(i,j,k)))
               enddo
            enddo
         enddo
         
         if (verbose == 1) then
            ntries = ntries+1
            print*,'-----------Retry', ntries
            print*,'residual =', rhsmax
            print*,'--------------------------------------------------'
         endif
         
         if (rhsmax.lt.rbh_tol) goto 110
c     
c     Compute stencil coefficients:
c     
         do k = 2,np-1
            do j = 2,nq-1
               do i = 2,ne-1
                  ac(i,j,k) = -2./deta**2-2./dq**2-2./(dphi**2*
     $                 sin(qgrd(j))**2)-0.25-0.875*Ksq(i,j,k)/(psisph(i,j,k)+
     $                 psip(i,j,k))**8
c     
                  ae(i,j,k) = 1./deta**2
c     
                  aw(i,j,k) = 1./deta**2
c         
                  an(i,j,k) = 1./dq**2 + 0.5/(tan(qgrd(j))*dq)
c     
                  as(i,j,k) = 1./dq**2 - 0.5/(tan(qgrd(j))*dq)
c     
                  aq(i,j,k) = 1./(sin(qgrd(j))**2*dphi**2)
c     
                  ab(i,j,k) = 1./(sin(qgrd(j))**2*dphi**2)
               enddo
            enddo
         enddo
c     
c     Apply boundary conditions:
c     
c     i=2:
c     
         do k = 3,np-2
            do j = 3,nq-2
               ae(2,j,k) = ae(2,j,k) + aw(2,j,k)
               aw(2,j,k) = 0.
c     
c     i=ne-1:
c     
               ac(ne-1,j,k) = ac(ne-1,j,k) + 4.*
     &              ae(ne-1,j,k)/(3.+deta)
               aw(ne-1,j,k) = aw(ne-1,j,k) -
     &              ae(ne-1,j,k)/(3.+deta)
               ae(ne-1,j,k) = 0.
            enddo
c     
c     j=2:
c     
            do i = 3,ne-2
               ac(i,2,k) = ac(i,2,k) + as(i,2,k)
               as(i,2,k) = 0.
c     
c     j=nq-1:
c     
               ac(i,nq-1,k) = ac(i,nq-1,k) +
     &              an(i,nq-1,k)
               an(i,nq-1,k) = 0.
            enddo
         enddo
c     
c     k=2:
c
         do j = 3,nq-2
            do i = 3,ne-2
               ac(i,j,2) = ac(i,j,2) + ab(i,j,2)
               ab(i,j,2) = 0.
c
c     k=np-1:
c
               ac(i,j,np-1) = ac(i,j,np-1) +
     &              aq(i,j,np-1)
               aq(i,j,np-1) = 0.
            enddo
         enddo
c     
c     i=2, j=2:
c
         do k = 3,np-2
            ae(2,2,k) = ae(2,2,k) + aw(2,2,k)
            ac(2,2,k) = ac(2,2,k) + as(2,2,k)
            aw(2,2,k) = 0.
            as(2,2,k) = 0.
c     
c     i=2, j=nq-1:
c
            ae(2,nq-1,k) = ae(2,nq-1,k) + aw(2,nq-1,k)
            ac(2,nq-1,k) = ac(2,nq-1,k) + an(2,nq-1,k)
            aw(2,nq-1,k) = 0.
            an(2,nq-1,k) = 0.
c     
c     i=ne-1, j=2:
c     
            ac(ne-1,2,k) = ac(ne-1,2,k) + 4.*ae(ne-1,2,k)/(3.+deta) +
     &           as(ne-1,2,k)
            aw(ne-1,2,k) = aw(ne-1,2,k) - ae(ne-1,2,k)/(3.+deta)
            ae(ne-1,2,k) = 0.
            as(ne-1,2,k) = 0.
c     
c     i=ne-1, j=nq-1:
c
            ac(ne-1,nq-1,k) = ac(ne-1,nq-1,k) + 4.*
     &           ae(ne-1,nq-1,k)/(3.+deta) + an(ne-1,nq-1,k)
            aw(ne-1,nq-1,k) = aw(ne-1,nq-1,k) -
     &           ae(ne-1,nq-1,k)/(3.+deta)
            ae(ne-1,nq-1,k) = 0.
            an(ne-1,nq-1,k) = 0.
         enddo
c     
c     i=2, k=2:
c     
         do j = 3,nq-2
            ae(2,j,2) = ae(2,j,2) + aw(2,j,2)
            ac(2,j,2) = ac(2,j,2) + ab(2,j,2)
            aw(2,j,2) = 0.
            ab(2,j,2) = 0.
c     
c     i=2, k=np-1:
c     
            ae(2,j,np-1) = ae(2,j,np-1) + aw(2,j,np-1)
            ac(2,j,np-1) = ac(2,j,np-1) + aq(2,j,np-1)
            aw(2,j,np-1) = 0.
            aq(2,j,np-1) = 0.
c
c     i=ne-1, k=2:
c     
            ac(ne-1,j,2) = ac(ne-1,j,2) + 4.*ae(ne-1,j,2)/(3.+deta) +
     &           ab(ne-1,j,2)
            aw(ne-1,j,2) = aw(ne-1,j,2) - ae(ne-1,j,2)/(3.+deta)
            ae(ne-1,j,2) = 0.
            ab(ne-1,j,2) = 0.
c     
c     i=ne-1, k=np-1:
c
            ac(ne-1,j,np-1) = ac(ne-1,j,np-1) + 4.*
     &           ae(ne-1,j,np-1)/(3.+deta) + aq(ne-1,j,np-1)
            aw(ne-1,j,np-1) = aw(ne-1,j,np-1) -
     &           ae(ne-1,j,np-1)/(3.+deta)
            ae(ne-1,j,np-1) = 0.
            aq(ne-1,j,np-1) = 0.
         enddo
c     
c     j=2, k=2:
c
         do i = 3,ne-2
            ac(i,2,2) = ac(i,2,2) + as(i,2,2) + ab(i,2,2)
            as(i,2,2) = 0.
            ab(i,2,2) = 0.
c
c     j=2, k=np-1:
c     
            ac(i,2,np-1) = ac(i,2,np-1) + as(i,2,np-1) +
     &           aq(i,2,np-1)
            as(i,2,np-1) = 0.
            aq(i,2,np-1) = 0.
c     
c  j=nq-1, k=2:
c     
            ac(i,nq-1,2) = ac(i,nq-1,2) + an(i,nq-1,2) +
     &           ab(i,nq-1,2)
            an(i,nq-1,2) = 0.
            ab(i,nq-1,2) = 0.
c
c  j=nq-1, k=np-1:
c
            ac(i,nq-1,np-1) = ac(i,nq-1,np-1) + an(i,nq-1,np-1) +
     &           aq(i,nq-1,np-1)
            an(i,nq-1,np-1) = 0.
            aq(i,nq-1,np-1) = 0.
         enddo
c
c  i=2, j=2, k=2:
c
         ae(2,2,2) = ae(2,2,2) + aw(2,2,2)
         ac(2,2,2) = ac(2,2,2) + as(2,2,2) + ab(2,2,2)
         aw(2,2,2) = 0.
         as(2,2,2) = 0.
         ab(2,2,2) = 0.
c
c  i=2, j=2, k=np-1:
c
         ae(2,2,np-1) = ae(2,2,np-1) + aw(2,2,np-1)
         ac(2,2,np-1) = ac(2,2,np-1) + as(2,2,np-1) + aq(2,2,np-1)
         aw(2,2,np-1) = 0.
         as(2,2,np-1) = 0.
         aq(2,2,np-1) = 0.
c
c  i=2, j=nq-1, k=2:
c
         ae(2,nq-1,2) = ae(2,nq-1,2) + aw(2,nq-1,2)
         ac(2,nq-1,2) = ac(2,nq-1,2) + an(2,nq-1,2) + ab(2,nq-1,2)
         aw(2,nq-1,2) = 0.
         an(2,nq-1,2) = 0.
         ab(2,nq-1,2) = 0.
c
c  i=2, j=nq-1, k=np-1:
c
         ae(2,nq-1,np-1) = ae(2,nq-1,np-1) + aw(2,nq-1,np-1)
         ac(2,nq-1,np-1) = ac(2,nq-1,np-1) + an(2,nq-1,np-1) + aq(2,nq-1,np-1)
         aw(2,nq-1,np-1) = 0.
         an(2,nq-1,np-1) = 0.
         aq(2,nq-1,np-1) = 0.
c
c  i=ne-1, j=2, k=2:
c
         ac(ne-1,2,2) = ac(ne-1,2,2) + 4.*ae(ne-1,2,2)/(3.+deta) +
     &        as(ne-1,2,2) + ab(ne-1,2,2)
         aw(ne-1,2,2) = aw(ne-1,2,2) - ae(ne-1,2,2)/(3.+deta)
         ae(ne-1,2,2) = 0.
         as(ne-1,2,2) = 0.
         ab(ne-1,2,2) = 0.
c
c  i=ne-1, j=2, k=np-1:
         ac(ne-1,2,np-1) = ac(ne-1,2,np-1) + 4.*ae(ne-1,2,np-1)/(3.+deta) +
     &        as(ne-1,2,np-1) + aq(ne-1,2,np-1)
         aw(ne-1,2,np-1) = aw(ne-1,2,np-1) - ae(ne-1,2,np-1)/(3.+deta)
         ae(ne-1,2,np-1) = 0.
         as(ne-1,2,np-1) = 0.
         aq(ne-1,2,np-1) = 0.
c
c  i=ne-1, j=nq-1, k=2:
c
         ac(ne-1,nq-1,2) = ac(ne-1,nq-1,2) + 4.*ae(ne-1,nq-1,2)/(3.+deta) +
     &        an(ne-1,nq-1,2) + ab(ne-1,nq-1,2)
         aw(ne-1,nq-1,2) = aw(ne-1,nq-1,2) - ae(ne-1,nq-1,2)/(3.+deta)
         ae(ne-1,nq-1,2) = 0.
         an(ne-1,nq-1,2) = 0.
         ab(ne-1,nq-1,2) = 0.
c
c  i=ne-1, j=nq-1, k=np-1:
c
         ac(ne-1,nq-1,np-1) = ac(ne-1,nq-1,np-1) +
     &        4.*ae(ne-1,nq-1,np-1)/(3.+deta) + an(ne-1,nq-1,np-1) +
     &        aq(ne-1,nq-1,np-1)
         aw(ne-1,nq-1,np-1) = aw(ne-1,nq-1,np-1) - ae(ne-1,nq-1,np-1)/
     $        (3.+deta)
         ae(ne-1,nq-1,np-1) = 0.
         an(ne-1,nq-1,np-1) = 0.
         aq(ne-1,nq-1,np-1) = 0.
c     
c     Call stab routine:
c     
         call rbicgst3d(ac,ae,aw,an,as,aq,ab,psipp,rhs,rbh_eps,rmax,ier,
     $        ne,nq,np)
         if (rmax.gt.1.0e-10) then
            write(*,*) '***WARNING: bicgst3d did not converge.'
         endif
         if (ier.eq.-1) then
            write(*,*) '***WARNING: ier=-1'
         endif
c     
c     Now, apply boundary conditions to psipp:
c     
         do k = 1,np
            do j = 1,nq
               psipp(1,j,k) = psipp(3,j,k)
               psipp(ne,j,k) = (4.*psipp(ne-1,j,k)-psipp(ne-2,j,k))/(3.+
     $              deta)
            enddo
            do i = 1,ne
               psipp(i,1,k) = psipp(i,2,k)
               psipp(i,nq,k) = psipp(i,nq-1,k)
            enddo
         enddo
         do j = 1,nq
            do i = 1,ne
               psipp(i,j,1) = psipp(i,j,2)
               psipp(i,j,np) = psipp(i,j,np-1)
            enddo
         enddo
c     
c     Update psip:
c     
         do k = 1,np
            do j = 1,nq
               do i = 1,ne
                  psip(i,j,k) = psip(i,j,k) + psipp(i,j,k)
               enddo
            enddo
         enddo
c     
      enddo
c     
      write(*,*) 'It did not converge.'
      stop
c     
 110  continue
      
      write(*,*) '--------------------------------------------------'
      print*,'Converge at Residual = ',rhsmax

c     
c     Here, compute the derivatives of the spherical conformal factor
c     
      do k = 1,np
         do j = 1,nq
            do i = 2,ne-1
               detapsisph(i,j,k)=0.5*(psip(i+1,j,k)-psip(i-1,j,k))/
     $              deta + sinh(0.5*etagrd(i))
            enddo
            detapsisph(1,j,k) = -detapsisph(3,j,k)
         enddo
      enddo
c     
      do k = 1,np
         do j = 2,nq-1
            do i = 1,ne
               dqpsisph(i,j,k)=0.5*(psip(i,j+1,k)-psip(i,j-1,k))/dq
            enddo
         enddo
         do i = 1,ne
            dqpsisph(i,1,k) = -dqpsisph(i,2,k)
            dqpsisph(i,nq,k) = -dqpsisph(i,nq-1,k)
         enddo
      enddo
c     
      do k = 2,np-1
         do j = 1,nq
            do i = 1,ne
               dphipsisph(i,j,k)=0.5*(psip(i,j,k+1)-psip(i,j,k-1))/
     $              dphi
            enddo
         enddo
      enddo
      do j = 1,nq
         do i = 1,ne
            dphipsisph(i,j,1) = -dphipsisph(i,j,2)
            dphipsisph(i,j,np) = -dphipsisph(i,j,np-1)
         enddo
      enddo
c     
      do k = 1,np
         do j = 1,nq
            do i = 2,ne-1
               detaetapsisph(i,j,k)=(psip(i+1,j,k)-2.*psip(i,j,k)+
     &              psip(i-1,j,k))/deta**2 + sqrt(0.25)*
     &              cosh(0.5*etagrd(i))
            enddo
            detaetapsisph(1,j,k) = detaetapsisph(3,j,k)
         enddo
      enddo
c     
      do k = 1,np
         do j = 2,nq-1
            do i = 1,ne
               detaqpsisph(i,j,k)=0.5*(detapsisph(i,j+1,k)-
     $              detapsisph(i,j-1,k))/dq
            enddo
         enddo
         do i = 1,ne
            detaqpsisph(i,1,k) = -detaqpsisph(i,2,k)
            detaqpsisph(i,nq,k) = -detaqpsisph(i,nq-1,k)
         enddo
      enddo
c     
      do k = 2,np-1
         do j = 1,nq
            do i = 1,ne
               detaphipsisph(i,j,k)=0.5*(detapsisph(i,j,k+1)-
     $              detapsisph(i,j,k-1))/dphi
            enddo
         enddo
      enddo
      do j = 1,nq
         do i = 1,ne
            detaphipsisph(i,j,1) = -detaphipsisph(i,j,2)
            detaphipsisph(i,j,np) = -detaphipsisph(i,j,np-1)
         enddo
      enddo
c     
      do k = 1,np
         do j = 2,nq-1
            do i = 1,ne
               dqqpsisph(i,j,k)=0.5*(dqpsisph(i,j+1,k)-
     $              dqpsisph(i,j-1,k))/dq
            enddo
         enddo
         do i = 1,ne
            dqqpsisph(i,1,k) = dqqpsisph(i,2,k)
            dqqpsisph(i,nq,k) = dqqpsisph(i,nq-1,k)
         enddo
      enddo
c     
      do k = 2,np-1
         do j = 1,nq
            do i = 1,ne
               dqphipsisph(i,j,k)=0.5*(dqpsisph(i,j,k+1)-
     $              dqpsisph(i,j,k-1))/dphi
            enddo
         enddo
      enddo
      do j = 1,nq
         do i = 1,ne
            dqphipsisph(i,j,1) = -dqphipsisph(i,j,2)
            dqphipsisph(i,j,np) = -dqphipsisph(i,j,np-1)
         enddo
      enddo
c     
      do k = 2,np-1
         do j = 1,nq
            do i = 1,ne
               dphiphipsisph(i,j,k)=0.5*(dphipsisph(i,j,k+1)-
     $              dphipsisph(i,j,k-1))/dphi
            enddo
         enddo
      enddo
      do j = 1,nq
         do i = 1,ne
            dphiphipsisph(i,j,1) = dphiphipsisph(i,j,2)
            dphiphipsisph(i,j,np) = dphiphipsisph(i,j,np-1)
         enddo
      enddo
c     
      do k = 1,np
         do j = 1,nq
            psisph(:,j,k)=psip(:,j,k)+2.0*cosh(0.5*etagrd)
         enddo
      enddo
c     
c     Now compute on the Cartesian coordinate.
c     
c     Compute eta,q,phi at the each points of cartesian grid
c     

      eta = 0.5d0 * log (x**2 + y**2 + z**2)
      abseta = abs (eta)
      q = atan2 (sqrt (x**2 + y**2), z)
      phi = atan2 (y, x)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx

               if(eta(i,j,k) .lt. 0)then
                  sign_eta(i,j,k) = -1.0d0
               else
                  sign_eta(i,j,k) = 1.0d0
               endif

            enddo
         enddo
      enddo
      
      npoints = nx*ny*nz

!     Parameter table and interpolator handles.
      options_string = "order = " // char(ichar('0') + interpolation_order)
      call Util_TableCreateFromString (param_table_handle, options_string)
      if (param_table_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot create parameter table for interpolator")
      endif

      call CCTK_InterpHandle (interp_handle, "uniform cartesian")
      if (interp_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
      endif

!     fill in the input/output arrays for the interpolator
      coord_origin(1) = etagrd(1)
      coord_origin(2) = qgrd(1)
      coord_origin(3) = phigrd(1)-pi
      coord_delta(1)  = etagrd(2) - etagrd(1)
      coord_delta(2)  = qgrd(2) - qgrd(1)
      coord_delta(3)  = phigrd(2) - phigrd(1)

      interp_coords(1) = CCTK_PointerTo(abseta)
      interp_coords(2) = CCTK_PointerTo(q)
      interp_coords(3) = CCTK_PointerTo(phi)

      in_array_dims(1) = ne
      in_array_dims(2) = nq
      in_array_dims(3) = np

      in_arrays(1) = CCTK_PointerTo(psisph)
      in_arrays(2) = CCTK_PointerTo(detapsisph)
      in_arrays(3) = CCTK_PointerTo(dqpsisph)
      in_arrays(4) = CCTK_PointerTo(dphipsisph)
      in_arrays(5) = CCTK_PointerTo(detaetapsisph)
      in_arrays(6) = CCTK_PointerTo(detaqpsisph)
      in_arrays(7) = CCTK_PointerTo(detaphipsisph)
      in_arrays(8) = CCTK_PointerTo(dqqpsisph)
      in_arrays(9) = CCTK_PointerTo(dqphipsisph)
      in_arrays(10) = CCTK_PointerTo(dphiphipsisph)

      out_arrays(1) = CCTK_PointerTo(psi3d)
      out_arrays(2) = CCTK_PointerTo(detapsi3d)
      out_arrays(3) = CCTK_PointerTo(dqpsi3d)
      out_arrays(4) = CCTK_PointerTo(dphipsi3d)
      out_arrays(5) = CCTK_PointerTo(detaetapsi3d)
      out_arrays(6) = CCTK_PointerTo(detaqpsi3d)
      out_arrays(7) = CCTK_PointerTo(detaphipsi3d)
      out_arrays(8) = CCTK_PointerTo(dqqpsi3d)
      out_arrays(9) = CCTK_PointerTo(dqphipsi3d)
      out_arrays(10) = CCTK_PointerTo(dphiphipsi3d)

      call CCTK_InterpLocalUniform (ierror, 3,
     $                              interp_handle, param_table_handle,
     $                              coord_origin, coord_delta,
     $                              npoints, type_codes(1), interp_coords,
     $                              10, in_array_dims, type_codes, in_arrays,
     $                              10, type_codes, out_arrays)

      psi = psi3d*exp(-0.5*eta)
      detapsi3d = sign_eta*detapsi3d
      detaqpsi3d = sign_eta*detaqpsi3d
      detaphipsi3d = sign_eta*detaphipsi3d

      do k=1,nz
         do j=1,ny
            do i=1,nx
c     psix = \partial psi / \partial x / psi
#include "psi_1st_deriv.x"
               
c     psixx = \partial^2\psi / \partial x^2 / psi
#include "psi_2nd_deriv.x"
            enddo
         enddo
      enddo
      
c     Conformal metric (here, we define conformally flat)
c     gxx = 1...
      
c     Derivatives of the metric
c     dxgxx = 1/2 \partial gxx / \partial x = 0...
c     
      do k= 1,nz
         do j= 1,ny
            do i= 1,nx
#include "gij.x"
            enddo
         enddo
      enddo
c     Extrinsic Curvture:Cactus only needs physical extrinsic curvature
      do k= 1,nz
         do j= 1,ny
            do i= 1,nx 
#include "cgauss.x"
#include "kij_odd.x"
            enddo
         enddo
      enddo
      
      kxx = kxx/psi**2
      kxy = kxy/psi**2
      kxz = kxz/psi**2
      kyy = kyy/psi**2
      kyz = kyz/psi**2
      kzz = kzz/psi**2

c     Set ADM mass
      i = ne-15
      adm = 0.0
      do k=2,np-1
         do j=2,nq-1
            adm=adm+(psisph(i,j,k)-(psisph(i+1,j,k)-psisph(i-1,j,k))/
     $           deta)*exp(0.5*etagrd(i))
         enddo
      enddo
      adm=adm/(nq-2)/(np-2)
      print *,'ADM mass: ',adm
      if (CCTK_Equals(initial_lapse,"schwarz")==1) then 
         write (*,*)"Initial with schwarzschild-like lapse"
         write (*,*)"using alp = (2.*r - adm)/(2.*r+adm)."
         alp = (2.*r - adm)/(2.*r+adm)
      endif
c
      print*,'Angular momentum parameter: a/m = ', byJ/adm**2
c
c     kerr functions
c     
c     isotropic coordinate for schwarzschild
c

      deallocate(ac,ae,aw,an,as,aq,ab,rhs,Ksq,psisph,psip,psipp,
     $     detapsisph,dqpsisph,dphipsisph,detaetapsisph,detaqpsisph,
     $     detaphipsisph,dqqpsisph,dqphipsisph,dphiphipsisph)
      deallocate(etagrd,qgrd,phigrd)
      
      return
      end

      
      
