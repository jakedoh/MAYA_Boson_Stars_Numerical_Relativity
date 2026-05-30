C $Header$

#include "cctk.h"
#include "cctk_Arguments.h"

      subroutine ComparisonMetric(CCTK_ARGUMENTS, 
     $     gxx_ex, gxy_ex, gxz_ex,
     $     gyy_ex, gyz_ex, gzz_ex)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      CCTK_REAL gxx_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL gxy_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL gxz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL gyy_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL gyz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL gzz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      
      integer i,j,k

C     Dummy arguments of subroutine boostrotdata: these are calculated
C     (at a point) but not used.

      CCTK_REAL hxxjunk, hyyjunk, hzzjunk, 
     $     hxyjunk, hyzjunk, hxzjunk,
     $     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $     dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $     dygxxjunk, dygyyjunk, dygzzjunk, 
     $     dygxyjunk, dygyzjunk, dygxzjunk,
     $     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $     dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $     alpjunk, axjunk, ayjunk, azjunk, 
     $     betaxjunk, betayjunk, betazjunk,
     $     bxxjunk, bxyjunk, bxzjunk, 
     $     byxjunk, byyjunk, byzjunk, 
     $     bzxjunk, bzyjunk, bzzjunk

C     Call boostrotdata pointwise. Most of what it calculates is
C     thrown away, variables ending in ...junk.
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
     $              gxx_ex(i,j,k), gyy_ex(i,j,k), gzz_ex(i,j,k), 
     $              gxy_ex(i,j,k), gyz_ex(i,j,k), gxz_ex(i,j,k),
     $              hxxjunk, hyyjunk, hzzjunk, 
     $              hxyjunk, hyzjunk, hxzjunk,
     $              dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $              dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $              dygxxjunk, dygyyjunk, dygzzjunk, 
     $              dygxyjunk, dygyzjunk, dygxzjunk,
     $              dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $              dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $              alpjunk, axjunk, ayjunk, azjunk, 
     $              betaxjunk, betayjunk, betazjunk,
     $              bxxjunk, bxyjunk, bxzjunk, 
     $              byxjunk, byyjunk, byzjunk, 
     $              bzxjunk, bzyjunk, bzzjunk)
            end do
         end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ComparisonCurvature(CCTK_ARGUMENTS,
     $     hxx_ex, hxy_ex, hxz_ex,
     $     hyy_ex, hyz_ex, hzz_ex)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      CCTK_REAL hxx_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL hxy_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL hxz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL hyy_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL hyz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL hzz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))

      integer i,j,k

      CCTK_REAL gxxjunk, gyyjunk, gzzjunk, 
     $     gxyjunk, gyzjunk, gxzjunk,
     $     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $     dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $     dygxxjunk, dygyyjunk, dygzzjunk, 
     $     dygxyjunk, dygyzjunk, dygxzjunk,
     $     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $     dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $     alpjunk, axjunk, ayjunk, azjunk, 
     $     betaxjunk, betayjunk, betazjunk,
     $     bxxjunk, bxyjunk, bxzjunk, 
     $     byxjunk, byyjunk, byzjunk, 
     $     bzxjunk, bzyjunk, bzzjunk

      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
     $              gxxjunk, gyyjunk, gzzjunk, 
     $              gxyjunk, gyzjunk, gxzjunk,
     $              hxx_ex(i,j,k), hyy_ex(i,j,k), hzz_ex(i,j,k), 
     $              hxy_ex(i,j,k), hyz_ex(i,j,k), hxz_ex(i,j,k),
     $              dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $              dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $              dygxxjunk, dygyyjunk, dygzzjunk, 
     $              dygxyjunk, dygyzjunk, dygxzjunk,
     $              dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $              dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $              alpjunk, axjunk, ayjunk, azjunk, 
     $              betaxjunk, betayjunk, betazjunk,
     $              bxxjunk, bxyjunk, bxzjunk, 
     $              byxjunk, byyjunk, byzjunk, 
     $              bzxjunk, bzyjunk, bzzjunk)
            end do
         end do
      end do
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Note the exact shift comes in even if the shift is not
c     allocated (eg if shift is "none"). In this case just
c     don't use it, since it won't be compared against acctk_lsh(2)thing.

      subroutine ComparisonGauge(CCTK_ARGUMENTS,
     $     alp_ex, betax_ex, betay_ex, betaz_ex)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      CCTK_REAL alp_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL betax_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL betay_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
      CCTK_REAL betaz_ex(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))

      integer i,j,k

      CCTK_REAL gxxjunk, gyyjunk, gzzjunk, 
     $     gxyjunk, gyzjunk, gxzjunk,
     $     dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $     dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $     dygxxjunk, dygyyjunk, dygzzjunk, 
     $     dygxyjunk, dygyzjunk, dygxzjunk,
     $     dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $     dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $     hxxjunk, hyyjunk, hzzjunk, 
     $     hxyjunk, hyzjunk, hxzjunk,
     $     axjunk, ayjunk, azjunk, 
     $     bxxjunk, bxyjunk, bxzjunk, 
     $     byxjunk, byyjunk, byzjunk, 
     $     bzxjunk, bzyjunk, bzzjunk

      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               call exactdata(x(i,j,k), y(i,j,k), z(i,j,k), cctk_time,
     $              gxxjunk, gyyjunk, gzzjunk, 
     $              gxyjunk, gyzjunk, gxzjunk,
     $              hxxjunk, hyyjunk, hzzjunk, 
     $              hxyjunk, hyzjunk, hxzjunk,
     $              dxgxxjunk, dxgyyjunk, dxgzzjunk, 
     $              dxgxyjunk, dxgyzjunk, dxgxzjunk,
     $              dygxxjunk, dygyyjunk, dygzzjunk, 
     $              dygxyjunk, dygyzjunk, dygxzjunk,
     $              dzgxxjunk, dzgyyjunk, dzgzzjunk, 
     $              dzgxyjunk, dzgyzjunk, dzgxzjunk,
     $              alp_ex(i,j,k), axjunk, ayjunk, azjunk, 
     $              betax_ex(i,j,k), betay_ex(i,j,k), betaz_ex(i,j,k),
     $              bxxjunk, bxyjunk, bxzjunk, 
     $              byxjunk, byyjunk, byzjunk, 
     $              bzxjunk, bzyjunk, bzzjunk)
            end do
         end do
      end do
      
      return
      end
