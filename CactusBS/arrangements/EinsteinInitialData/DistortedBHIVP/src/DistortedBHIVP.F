c     /*@@
c     @file      DistortedBHIVP.F
c     @date      
c     @author    
c     @desc 
c     
c     @enddesc 
c@@   */
      
#include "cctk.h" 

#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      
c     /*@@
c     @routine    DistortedBHIVP
c     @date       
c     @author     
c     @desc 
c     
c     @enddesc 
c     @calls     
c     @calledby   
c     @history 
c     
c     @endhistory 
c@@   */

      subroutine DistortedBHIVP(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      real*8 :: deta,dq,dphi
      real*8, allocatable :: ac(:,:,:),ae(:,:,:),aw(:,:,:),an(:,:,:),
     $     as(:,:,:),aq(:,:,:),ab(:,:,:),rhs(:,:,:),
     $     qf(:,:,:),qfetaeta(:,:,:),qfqq(:,:,:),qfphi(:,:,:),
     $     qfphiphi(:,:,:),
     $     psisph(:,:,:),psiprim(:,:,:),detapsisph(:,:,:),
     $     dqpsisph(:,:,:),dphipsisph(:,:,:),detaetapsisph(:,:,:),
     $     detaqpsisph(:,:,:),detaphipsisph(:,:,:),dqqpsisph(:,:,:),
     $     dqphipsisph(:,:,:),dphiphipsisph(:,:,:)
      real*8, allocatable :: etagrd(:),qgrd(:),phigrd(:)
      real*8 o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,o17,
     $     o18,o19,o20,o21,o22,o23,o24,o25,o26,o27,o28,o29,o30,o31,o32,
     $     o33,o34,o35,o36,o37,o38,o39,o40,o41,o42,o43,o44,o45,o46,o47,
     $     o48,o49,o50,o51,o52,o53,o54,o55,o56,o57,o58,o59,o60,o61,o62,
     $     o63,o64
      real*8 rmax,adm
      real*8 dbh_eps
      real*8 pi
      integer :: ne, nq, np
      integer :: nx,ny,nz
      integer i,j,k,ier,nquads,nocts
      integer npoints,ierror
      integer make_conformal_derivs

      integer       param_table_handle, interp_handle
      character(30) options_string
      CCTK_REAL,    dimension(3)  :: coord_origin, coord_delta
      CCTK_POINTER, dimension(3)  :: interp_coords
      CCTK_INT,     dimension(3 ) :: in_array_dims
      CCTK_POINTER, dimension(10) :: in_arrays, out_arrays
      CCTK_INT,     dimension(10), parameter :: type_codes = CCTK_VARIABLE_REAL

      dbh_eps = tolerance

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

c     DONT use integer*4
      ne = neta
      nq = ntheta
      np = nphi

c     Set up the grid spacings
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)      
      
c     Distorted Schwarzchild BH parameters
      
      print *,"Brill wave + Distorted BH solve"
      write(*,123)amp,eta0,c,sigma,n
      print*,'etamax=',etamax
 123  format(1x, 'Pars: amp',f8.5,' eta0',f8.5,' c',f8.5,' sigma',f8.5,' n',i3)
      
c     Sovle on this sized cartesian grid
c     3D grid size NE x NT x NP
c     Add 2 zones for eta coordinate and 4 for theta
c     and phi coordenate.
      ne = ne + 2
      nq = nq + 2
      np = np + 2
c     
      allocate(ac(ne,nq,np),ae(ne,nq,np),aw(ne,nq,np),an(ne,nq,np),
     $     as(ne,nq,np),aq(ne,nq,np),ab(ne,nq,np),rhs(ne,nq,np),
     $     qf(ne,nq,np),qfetaeta(ne,nq,np),qfqq(ne,nq,np),
     $     qfphi(ne,nq,np),qfphiphi(ne,nq,np),
     $     psisph(ne,nq,np),psiprim(ne,nq,np),detapsisph(ne,nq,np),
     $     dqpsisph(ne,nq,np),dphipsisph(ne,nq,np),
     $     detaetapsisph(ne,nq,np),detaqpsisph(ne,nq,np),
     $     detaphipsisph(ne,nq,np),dqqpsisph(ne,nq,np),
     $     dqphipsisph(ne,nq,np),dphiphipsisph(ne,nq,np))
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
c Initialize q-function and its derivatives: should be generalized
c     
      do k = 1,np
         do j = 1,nq
            do i = 1,ne
#include "qfunc.x"
            enddo
         enddo
      enddo
      
c     
c     Initialize psi to the Schwarzschild solution:
c     
      psiprim = 0.
      
      do k = 1,np
         do j = 1,nq
            do i = 1,ne
               psisph(i,j,k) = 2.*cosh(0.5*etagrd(i))
            enddo
         enddo
      enddo
c     
c     Initialize stencil coefficients:
c     
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
               ac(i,j,k) = -2./deta**2-2./dq**2-2.*exp(2.*
     &              qf(i,j,k))/(dphi**2*
     &              sin(qgrd(j))**2)+0.25*
     &              (qfetaeta(i,j,k)+qfqq(i,j,k)+2.*
     &              exp(2.*qf(i,j,k))*qfphiphi(i,j,k)/
     &              sin(qgrd(j))**2+3.*exp(2.*
     &              qf(i,j,k))*qfphi(i,j,k)**2/
     &              sin(qgrd(j))**2-1.)
               ae(i,j,k) = 1./deta**2
               aw(i,j,k) = 1./deta**2
               an(i,j,k) = 1./dq**2+0.5/(dq*tan(qgrd(j)))
               
               as(i,j,k) = 1./dq**2-0.5/(dq*tan(qgrd(j)))
               
               aq(i,j,k) = exp(2.*qf(i,j,k))/(dphi**2*
     &              sin(qgrd(j))**2)+exp(2.*
     &              qf(i,j,k))*qfphi(i,j,k)/(dphi*
     &              sin(qgrd(j))**2)
               
               ab(i,j,k) = exp(2.*qf(i,j,k))/(dphi**2*
     &              sin(qgrd(j))**2)-exp(2.*
     &              qf(i,j,k))*qfphi(i,j,k)/(dphi*
     &              sin(qgrd(j))**2)
               
               rhs(i,j,k) = -0.25*(qfetaeta(i,j,k)+
     &              qfqq(i,j,k)+2.*exp(2.*qf(i,j,k))*
     &              qfphiphi(i,j,k)/sin(qgrd(j))**2+3.*
     &              exp(2.*qf(i,j,k))*qfphi(i,j,k)**2/
     &              sin(qgrd(j))**2)*psisph(i,j,k)
            enddo
         enddo
      enddo
c
c Apply boundary conditions to the faces of the cube:
c
c    i=2:
      do k = 3,np-2
         do j = 3,nq-2
            ae(2,j,k) = ae(2,j,k) + aw(2,j,k)
            aw(2,j,k) = 0.
c
c    i=ne-1:
            ac(ne-1,j,k) = ac(ne-1,j,k)+4.*ae(ne-1,j,k)/(3.+deta)
            aw(ne-1,j,k) = aw(ne-1,j,k) - ae(ne-1,j,k)/(3.+deta)
            ae(ne-1,j,k) = 0.
         enddo
      enddo
c
c    j=2:
      do k = 3,np-2
         do i = 3,ne-2
            ac(i,2,k) = ac(i,2,k) + as(i,2,k)
            as(i,2,k) = 0.
c
c    j=nq-1:
            ac(i,nq-1,k) = ac(i,nq-1,k) + an(i,nq-1,k)
            an(i,nq-1,k) = 0.
         enddo
      enddo
c
c    k=2:
      do j = 3,nq-2
         do i = 3,ne-2
            ac(i,j,2) = ac(i,j,2) + ab(i,j,2)
            ab(i,j,2) = 0.
c
c    k=np-1:
            ac(i,j,np-1) = ac(i,j,np-1) + aq(i,j,np-1)
            aq(i,j,np-1) = 0.
         enddo
      enddo
c
c Apply boundary conditions to the edges of the cube:
c
c    i=2, j=2:
      do k = 3,np-2
         ae(2,2,k) = ae(2,2,k) + aw(2,2,k)
         ac(2,2,k) = ac(2,2,k) + as(2,2,k)
         aw(2,2,k) = 0.
         as(2,2,k) = 0.
c
c    i=ne-1, j=2:
         aw(ne-1,2,k) = aw(ne-1,2,k) - ae(ne-1,2,k)/(3.+deta)
         ac(ne-1,2,k) = ac(ne-1,2,k) + as(ne-1,2,k) +
     &        4.*ae(ne-1,2,k)/(3.+deta)
         ae(ne-1,2,k) = 0.
         as(ne-1,2,k) = 0.
c
c    i=2, j=nq-1:
         ae(2,nq-1,k) = ae(2,nq-1,k) + aw(2,nq-1,k)
         ac(2,nq-1,k) = ac(2,nq-1,k) + an(2,nq-1,k)
         aw(2,nq-1,k) = 0.
         an(2,nq-1,k) = 0.
c
c    i=ne-1, j=nq-1:
         aw(ne-1,nq-1,k) = aw(ne-1,nq-1,k) - ae(ne-1,nq-1,k)/
     &        (3.+deta)
         ac(ne-1,nq-1,k) = ac(ne-1,nq-1,k) + an(ne-1,nq-1,k) +
     &        4.*ae(ne-1,nq-1,k)/(3.+deta)
         ae(ne-1,nq-1,k) = 0.
         an(ne-1,nq-1,k) = 0.
      enddo
c
c    i=2, k=2:
      do j = 3,nq-2
         ae(2,j,2) = ae(2,j,2) + aw(2,j,2)
         ac(2,j,2) = ac(2,j,2) + ab(2,j,2)
         aw(2,j,2) = 0.
         ab(2,j,2) = 0.
c
c    i=ne-1, k=2:
         aw(ne-1,j,2) = aw(ne-1,j,2) - ae(ne-1,j,2)/(3.+deta)
         ac(ne-1,j,2) = ac(ne-1,j,2) + ab(ne-1,j,2) +
     &        4.*ae(ne-1,j,2)/(3.+deta)
         ae(ne-1,j,2) = 0.
         ab(ne-1,j,2) = 0.
c
c    i=2, k=np-1:
         ae(2,j,np-1) = ae(2,j,np-1) + aw(2,j,np-1)
         ac(2,j,np-1) = ac(2,j,np-1) + aq(2,j,np-1)
         aw(2,j,np-1) = 0.
         aq(2,j,np-1) = 0.
c
c    i=ne-1, k=np-1:
         aw(ne-1,j,np-1) = aw(ne-1,j,np-1) - ae(ne-1,j,np-1)/
     &        (3.+deta)
         ac(ne-1,j,np-1) = ac(ne-1,j,np-1) + aq(ne-1,j,np-1) +
     &        4.*ae(ne-1,j,np-1)/(3.+deta)
         ae(ne-1,j,np-1) = 0.
         aq(ne-1,j,np-1) = 0.
      enddo
c
c    j=2, k=2:
      do i = 3,ne-2
         ac(i,2,2) = ac(i,2,2) + as(i,2,2) + ab(i,2,2)
         as(i,2,2) = 0.
         ab(i,2,2) = 0.
c
c    j=nq-1, k=2:
         ac(i,nq-1,2) = ac(i,nq-1,2) + an(i,nq-1,2) +
     &        ab(i,nq-1,2)
         an(i,nq-1,2) = 0.
         ab(i,nq-1,2) = 0.
c
c    j=2, k=np-1:
         ac(i,2,np-1) = ac(i,2,np-1) + as(i,2,np-1) +
     &        aq(i,2,np-1)
         as(i,2,np-1) = 0.
         aq(i,2,np-1) = 0.
c
c    j=nq-1, k=np-1:
         ac(i,nq-1,np-1) = ac(i,nq-1,np-1) + an(i,nq-1,np-1) +
     &        aq(i,nq-1,np-1)
         an(i,nq-1,np-1) = 0.
         aq(i,nq-1,np-1) = 0.
      enddo
c
c Apply boundary conditions to the corners of the cube:
c
c    i=2, j=2, k=2:
      ae(2,2,2) = ae(2,2,2) + aw(2,2,2)
      ac(2,2,2) = ac(2,2,2) + as(2,2,2) + ab(2,2,2)
      aw(2,2,2) = 0.
      as(2,2,2) = 0.
      ab(2,2,2) = 0.
c
c    i=ne-1, j=2, k=2:
      aw(ne-1,2,2) = aw(ne-1,2,2) - ae(ne-1,2,2)/(3.+deta)
      ac(ne-1,2,2) = ac(ne-1,2,2) + as(ne-1,2,2) + ab(ne-1,2,2) +
     &     4.*ae(ne-1,2,2)/(3.+deta)
      ae(ne-1,2,2) = 0.
      as(ne-1,2,2) = 0.
      ab(ne-1,2,2) = 0.
c
c    i=2, j=nq-1, k=2:
      ae(2,nq-1,2) = ae(2,nq-1,2) + aw(2,nq-1,2)
      ac(2,nq-1,2) = ac(2,nq-1,2) + an(2,nq-1,2) + ab(2,nq-1,2)
      aw(2,nq-1,2) = 0.
      an(2,nq-1,2) = 0.
      ab(2,nq-1,2) = 0.
c
c    i=2, j=2, k=np-1:
      ae(2,2,np-1) = ae(2,2,np-1) + aw(2,2,np-1)
      ac(2,2,np-1) = ac(2,2,np-1) + as(2,2,np-1) + aq(2,2,np-1)
      aw(2,2,np-1) = 0.
      as(2,2,np-1) = 0.
      aq(2,2,np-1) = 0.
c
c    i=ne-1, j=nq-1, k=2:
      aw(ne-1,nq-1,2) = aw(ne-1,nq-1,2) - ae(ne-1,nq-1,2)/(3.+deta)
      ac(ne-1,nq-1,2) = ac(ne-1,nq-1,2) + an(ne-1,nq-1,2) + ab(ne-1,nq-1,2) +
     &     4.*ae(ne-1,nq-1,2)/(3.+deta)
      ae(ne-1,nq-1,2) = 0.
      an(ne-1,nq-1,2) = 0.
      ab(ne-1,nq-1,2) = 0.
c
c    i=ne-1, j=2, k=np-1:
      aw(ne-1,2,np-1) = aw(ne-1,2,np-1) - ae(ne-1,2,np-1)/(3.+deta)
      ac(ne-1,2,np-1) = ac(ne-1,2,np-1) + as(ne-1,2,np-1) + aq(ne-1,2,np-1) +
     &     4.*ae(ne-1,2,np-1)/(3.+deta)
      ae(ne-1,2,np-1) = 0.
      as(ne-1,2,np-1) = 0.
      aq(ne-1,2,np-1) = 0.
c
c    i=2, j=nq-1, k=np-1:
      ae(2,nq-1,np-1) = ae(2,nq-1,np-1) + aw(2,nq-1,np-1)
      ac(2,nq-1,np-1) = ac(2,nq-1,np-1) + an(2,nq-1,np-1) + aq(2,nq-1,np-1)
      aw(2,nq-1,np-1) = 0.
      an(2,nq-1,np-1) = 0.
      aq(2,nq-1,np-1) = 0.
c
c    i=ne-1, j=nq-1, k=np-1:
      aw(ne-1,nq-1,np-1) = aw(ne-1,nq-1,np-1) - ae(ne-1,nq-1,np-1)/(3.+deta)
      ac(ne-1,nq-1,np-1) = ac(ne-1,nq-1,np-1) + an(ne-1,nq-1,np-1) +
     &     aq(ne-1,nq-1,np-1) + 4.*ae(ne-1,nq-1,np-1)/(3.+deta)
      ae(ne-1,nq-1,np-1) = 0.
      an(ne-1,nq-1,np-1) = 0.
      aq(ne-1,nq-1,np-1) = 0.
c
c Solve for psi:
c
      call bicgst3d(ac,ae,aw,an,as,aq,ab,psiprim,rhs,dbh_eps,rmax,ier,
     $     ne,nq,np)
c
      if (rmax.gt.dbh_eps) then
         write(*,*) '***WARNING: bicgst3d did not converge.'
      endif
      if (ier.eq.-1) then
         write(*,*) '***WARNING: ier=-1'
      endif
      print *,'psiprim = ',maxval(psiprim),' ',minval(psiprim)
c
c Now, apply boundary conditions to psiprim:
c
      do k = 1,np
         do j = 1,nq
            psiprim(1,j,k) = psiprim(3,j,k)
            psiprim(ne,j,k) = (4.*psiprim(ne-1,j,k)-psiprim(ne-2,j,k))/
     $           (3.+deta)
         enddo
      enddo
      do k = 1,np
         do i = 1,ne
            psiprim(i,1,k) = psiprim(i,2,k)
            psiprim(i,nq,k) = psiprim(i,nq-1,k)
         enddo
      enddo
      do j = 1,nq
         do i = 1,ne
            psiprim(i,j,1) = psiprim(i,j,2)
            psiprim(i,j,np) = psiprim(i,j,np-1)
         enddo
      enddo
c     
c     Here, compute the derivatives of the spherical conformal factor
c
c      goto 110 

      do k = 1,np
         do j = 1,nq
            do i = 2,ne-1
               detapsisph(i,j,k)=0.5*(psiprim(i+1,j,k)-psiprim(i-1,j,k))
     $              /deta + sinh(0.5*etagrd(i))
            enddo
            detapsisph(1,j,k) = -detapsisph(3,j,k)
         enddo
      enddo
c     
      do k = 1,np
         do j = 2,nq-1
            do i = 1,ne
               dqpsisph(i,j,k)=0.5*(psiprim(i,j+1,k)-psiprim(i,j-1,k))/
     $              dq
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
               dphipsisph(i,j,k)=0.5*(psiprim(i,j,k+1)-psiprim(i,j,k-1))
     $              /dphi
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
               detaetapsisph(i,j,k)=(psiprim(i+1,j,k)-2.*psiprim(i,j,k)+
     &              psiprim(i-1,j,k))/deta**2+sqrt(0.25)*
     &              cosh(0.5*etagrd(i))
            enddo
            detaetapsisph(1,j,k) = detaetapsisph(3,j,k)
         enddo
      enddo
c     
      do k = 1,np
         do j = 2,nq-1
            do i = 1,ne
               detaqpsisph(i,j,k)=0.5*(detapsisph(i,j+1,k)-detapsisph(i,
     $              j-1,k))/dq
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
               detaphipsisph(i,j,k)=0.5*(detapsisph(i,j,k+1)-detapsisph(
     $              i,j,k-1))/dphi
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
               dqqpsisph(i,j,k)=0.5*(dqpsisph(i,j+1,k)-dqpsisph(i,j-1,k))/
     $              dq
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
               dqphipsisph(i,j,k)=0.5*(dqpsisph(i,j,k+1)-dqpsisph(i,j,k-1))/
     $              dphi
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
            psisph(:,j,k)=psiprim(:,j,k)+2.0*cosh(0.5*etagrd)
         enddo
      enddo
      
c     
c     Now compute on the Cartesian coordinate.
c     
c     Compute eta,q,phi at the each points of cartesian grid

      eta = 0.5d0*log(x**2+y**2+z**2)
      abseta = abs(eta)
      q = atan2(sqrt(x**2+y**2),z)
      phi = atan2(y,x)

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
      
c     Conformal metric
c     gxx = ...
      
c     Derivatives of the metric
c     dxgxx = 1/2 \partial gxx / \partial x
c     
      do k=1,nz
         do j=1,ny
            do i=1,nx
#include "gij.x"
            enddo
         enddo
      enddo
      
c     Courvature
      kxx = 0.0d0
      kxy = 0.0d0
      kxz = 0.0d0
      kyy = 0.0d0
      kyz = 0.0d0
      kzz = 0.0d0
      
c 110  continue
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
      
      deallocate(ac,ae,aw,an,as,aq,ab,rhs,
     $     qf,qfetaeta,qfqq,qfphi,qfphiphi,
     $     psisph,psiprim,detapsisph,dqpsisph,dphipsisph,detaetapsisph,
     $     detaqpsisph,detaphipsisph,dqqpsisph,dqphipsisph,
     $     dphiphipsisph)
      deallocate(etagrd,qgrd,phigrd)
      
      return
      end
