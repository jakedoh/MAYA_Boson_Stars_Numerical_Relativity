
#include "cctk.h"

      MODULE D3_to_D2_int

c     ------------------------------------------------------------------

      INTERFACE

      SUBROUTINE D3_to_D2(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,
     &        Psi_power,origin,myproc,interpolation_operator,interpolation_order,
     &        Dx,Dy,Dz,Psi,
     &        g00,gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        x,y,z,eta,Nt,Np,nx,ny,nz,theta,phi,Psis,g00s,gxxs,gxys,
     &        gxzs,gyys,gyzs,gzzs,dPsis,dgxxs,dgxys,dgxzs,dgyys,dgyzs,
     &        dgzzs,ADMmass_int1,ADMmass_int2,
     &        momentum_int1,momentum_int2,momentum_int3,
     &        spin_int1,spin_int2,spin_int3,Extract_temp3d)

      IMPLICIT NONE
      CCTK_POINTER :: cctkGH
      CCTK_INT,INTENT(IN) :: 
     &     conformal_state
      INTEGER,INTENT(IN) :: 
     &     myproc,Psi_power
      CCTK_INT, INTENT(IN) ::
     &     Nt,Np,nx,ny,nz,do_momentum,do_spin,interpolation_order
      CCTK_REAL,INTENT(IN) :: 
     &     origin(3),Dx,Dy,Dz,eta
      CCTK_REAL,INTENT(IN),DIMENSION(Nt) :: theta
      CCTK_REAL,INTENT(IN),DIMENSION(Np) :: phi
      CCTK_REAL,INTENT(IN),DIMENSION(nx) :: x
      CCTK_REAL,INTENT(IN),DIMENSION(ny) :: y
      CCTK_REAL,INTENT(IN),DIMENSION(nz) :: z
      CCTK_REAL,INTENT(IN),DIMENSION(nx,ny,nz) ::
     &     Psi,g00,gxx,gxy,gxz,gyy,gyz,gzz,
     &     hxx,hxy,hxz,hyy,hyz,hzz
      INTEGER,INTENT(IN) :: do_ADMmass(2)
      CCTK_STRING,INTENT(IN) ::
     &     interpolation_operator

      CCTK_REAL,INTENT(INOUT),DIMENSION(nx,ny,nz) ::
     &     Extract_temp3d

      CCTK_REAL,INTENT(OUT),DIMENSION(Nt,Np) :: 
     &     Psis,g00s,gxxs,gxys,gxzs,gyys,
     &     gyzs,gzzs,dPsis,dgxxs,dgxys,dgxzs,dgyys,dgyzs,dgzzs,
     &     ADMmass_int1,ADMmass_int2,
     &     momentum_int1,momentum_int2,momentum_int3,
     &     spin_int1,spin_int2,spin_int3
      END SUBROUTINE

      END INTERFACE

      END MODULE D3_to_D2_int
