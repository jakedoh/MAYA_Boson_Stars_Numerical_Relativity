c=======================================================================
c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c by Hisaaki Shinkai   shinkai@wurel.wustl.edu  19980603
c-----------------------------------------------------------------------
c  This is for maximally charged multi BH solutions such as 
c  Majumdar-Papapetrou (1947) or Kastor-Traschen (1993) solution.
c  See also doc/KTsol.tex for brief review of this solution. 
c-----------------------------------------------------------------------
c  For usage:  in your par file
c     Exact::exact_model = "multiBH"
c     Exact::multi_BH__Hubble       = 0.0   # 0.0 means MP solution
c     Exact::multi_BH__nBH          = 2     # number of BHs (upto 4, currently)
c     m_bh1,multi_BH__x1x,multi_BH__x1y,multi_BH__x1z =  1.0, -2.0,0.0,0.0  # masses and 
c     m_bh2,multi_BH__x2x,multi_BH__x2y,multi_BH__x2z =  1.0,  2.0,0.0,0.0  #   locations
c-----------------------------------------------------------------------
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__multi_BH(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local static variables
      logical firstcall
      CCTK_REAL kt_xbh(10),kt_ybh(10),kt_zbh(10),kt_mbh(10)
      data firstcall /.true./
      save firstcall,kt_xbh,kt_ybh,kt_zbh,kt_mbh

c local variables
      CCTK_REAL kt_r, kt_aa, kt_omega
      integer i

c constants
      CCTK_REAL zero,one
      parameter (zero=0.0d0, one=1.0d0)

C This model does not set the stress-energy tensor
C FIXME: should it?  I.e. isnt there a nonzero Maxwell tensor here?
      Tmunu_flag = .false.

c     Get parameters of the exact solution.

      if (firstcall) then

         write(*,*) ' welcome to Kastor-Traschen (Majumdar-Papapetrou)'

         if(multi_BH__nBH.ge.1) then
            kt_xbh(1) = multi_BH__x1
            kt_ybh(1) = multi_BH__y1
            kt_zbh(1) = multi_BH__z1
            kt_mbh(1) = multi_BH__mass1
         endif

         if(multi_BH__nBH.ge.2) then
            kt_xbh(2) = multi_BH__x2
            kt_ybh(2) = multi_BH__y2
            kt_zbh(2) = multi_BH__z2
            kt_mbh(1) = multi_BH__mass2
         endif

         if(multi_BH__nBH.ge.3) then
            kt_xbh(3) = multi_BH__x3
            kt_ybh(3) = multi_BH__y3
            kt_zbh(3) = multi_BH__z3
            kt_mbh(1) = multi_BH__mass3
         endif

         if(multi_BH__nBH.ge.4) then
            kt_xbh(4) = multi_BH__x4
            kt_ybh(4) = multi_BH__y4
            kt_zbh(4) = multi_BH__z4
            kt_mbh(1) = multi_BH__mass4
         endif

         write(*,*) 'time=',t
         write(*,*) '   mass    BH(x,y,z) '

         do i=1,multi_BH__nBH
            write(*,'(4e12.3)') kt_mbh(i),kt_xbh(i),kt_ybh(i),kt_zbh(i)
         enddo

         firstcall = .false.

      end if

      kt_aa=exp(multi_BH__Hubble*t)
      kt_omega=1.0

      do i=1,multi_BH__nBH
         kt_r=sqrt((x-kt_xbh(i))**2+(y-kt_ybh(i))**2+(z-kt_zbh(i))**2)
         kt_omega= kt_omega+kt_mbh(i)/kt_aa/kt_r
      enddo
      
c     write(*,*) kt_omega,kt_aa

      gdtt = -1.0/kt_omega**2
      gdtx = zero
      gdty = zero
      gdtz = zero
      gdxx = (kt_aa*kt_omega)**2
      gdyy = (kt_aa*kt_omega)**2
      gdzz = (kt_aa*kt_omega)**2
      gdxy = zero
      gdyz = zero
      gdzx = zero

      gutt = one/gdtt 
      gutx = zero
      guty = zero
      gutz = zero
      guxx = one/gdxx
      guyy = one/gdyy
      guzz = one/gdzz
      guxy = zero
      guyz = zero
      guzx = zero

      return
      end
