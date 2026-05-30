
#include "cctk.h"

c     ==================================================================

      SUBROUTINE cartesian_to_spherical(theta,phi,r,gxxs,gxys,gxzs,
     &     gyys,gyzs, gzzs,grrs,grts,grps,gtts,gtps,gpps,dgxxs,dgxys,
     &     dgxzs,dgyys,dgyzs,dgzzs,dgtts,dgtps,dgpps)

c     ------------------------------------------------------------------
c
c     Convert components of a symmetric 2nd rank 3D tensor on a 2-sphere 
c     from Cartesian to spherical polar coordinates
c
c                (x,y,z)        ---->         (r,t,p)
c
c           ( gxx  gxy gxz )             ( grr  grt grp )
c           (  .   gyy gyz )    ---->    (  .   gpp gtp )   
c           (  .    .  gzz )             (  .    .  gpp ) 
c
c     Also convert radial derivatives of the tensor components
c
c     ------------------------------------------------------------------

      IMPLICIT NONE

c     Input variables
      CCTK_REAL,INTENT(IN) ::  
     &    r
      CCTK_REAL,INTENT(IN),DIMENSION (:) ::  
     &    theta,phi
      CCTK_REAL,INTENT(IN),DIMENSION (:,:) :: 
     &    gxxs,gxys,gxzs,gyys,gyzs,gzzs,dgxxs,dgxys,dgxzs,dgyys,dgyzs,
     &    dgzzs

c     Output variables
      CCTK_REAL,INTENT(OUT),DIMENSION (:,:) :: 
     &    grrs,grts,grps,gtts,gtps,
     &    gpps,dgtts,dgtps,dgpps

c     Local variables
      INTEGER ::
     &    it,ip,Nt,Np
      CCTK_REAL :: 
     &    ct,st,cp,sp,ct2,st2,cp2,sp2,r2,gxx,gxy,gxz,gyy,gyz,gzz,dgxx,
     &    dgxy,dgxz,dgyy,dgyz,dgzz
      CCTK_REAL,PARAMETER ::
     &    two = 2.0D0

c     ------------------------------------------------------------------

      Nt = SIZE(theta)
      Np = SIZE(phi)

      r2 = r*r

      DO it = 1, Nt
 
        ct = COS(theta(it)) ; ct2 = ct*ct
        st = SIN(theta(it)) ; st2 = st*st

          DO ip = 1, Np

          cp = COS(phi(ip)) ; cp2 = cp*cp
          sp = SIN(phi(ip)) ; sp2 = sp*sp

          gxx = gxxs(it,ip) ; dgxx = dgxxs(it,ip)
          gxy = gxys(it,ip) ; dgxy = dgxys(it,ip)
          gxz = gxzs(it,ip) ; dgxz = dgxzs(it,ip)
          gyy = gyys(it,ip) ; dgyy = dgyys(it,ip)
          gyz = gyzs(it,ip) ; dgyz = dgyzs(it,ip)
          gzz = gzzs(it,ip) ; dgzz = dgzzs(it,ip)

          grrs(it,ip) = st2*cp2*gxx+st2*sp2*gyy+ct2*gzz
     &       +two*st2*cp*sp*gxy+two*st*cp*ct*gxz+two*st*ct*sp*gyz

          grts(it,ip) = r*(st*cp2*ct*gxx+two*st*ct*sp*cp*gxy
     &       +cp*(ct2-st2)*gxz+st*sp2*ct*gyy+sp*(ct2-st2)*gyz-ct*st*gzz)

          grps(it,ip) = r*st*(-st*sp*cp*gxx-st*(sp2-cp2)*gxy-sp*ct*gxz
     &       +st*sp*cp*gyy+ct*cp*gyz)
     
          gtts(it,ip) = r2*(ct2*cp2*gxx+two*ct2*sp*cp*gxy
     &       -two*st*ct*cp*gxz+ct2*sp2*gyy-two*st*sp*ct*gyz+st2*gzz)
     
          gtps(it,ip) = r2*st*(-cp*sp*ct*gxx-ct*(sp2-cp2)*gxy
     &       +st*sp*gxz+cp*sp*ct*gyy-st*cp*gyz)
     
          gpps(it,ip) = r2*st2*(sp2*gxx-two*cp*sp*gxy+cp2*gyy)
     
          dgtts(it,ip) = two/r*gtts(it,ip)+r2*(ct2*cp2*dgxx
     &       +two*ct2*sp*cp*dgxy-two*st*ct*cp*dgxz+ct2*sp2*dgyy
     &       -two*st*sp*ct*dgyz+st2*dgzz)

          dgtps(it,ip) = two/r*gtps(it,ip)+r2*st*(-cp*sp*ct*dgxx
     &       -ct*(sp2-cp2)*dgxy+st*sp*dgxz+cp*sp*ct*dgyy-st*cp*dgyz)
     
          dgpps(it,ip) = two/r*gpps(it,ip) +r2*st2*(sp2*dgxx
     &       -two*cp*sp*dgxy+cp2*dgyy)

        END DO
                                         
      END DO


      END SUBROUTINE cartesian_to_spherical










