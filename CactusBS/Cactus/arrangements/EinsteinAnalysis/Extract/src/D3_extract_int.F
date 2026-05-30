
#include "cctk.h"

      MODULE D3_extract_int

c     ------------------------------------------------------------------

      INTERFACE

      SUBROUTINE D3_extract(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,igrid,
     &        origin,myproc,interpolation_operator,interpolation_order,
     &        Nt,Np,nx,ny,nz,all_modes,l,m,x,y,z,Dx,Dy,Dz,Psi_power,Psi,
     &        g00,gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &        eta,ADMmass,momentum,spin,mass,rsch,Qodd,Qeven,
     &        Extract_temp3d,dtaudt)

      IMPLICIT NONE

      CCTK_POINTER :: cctkGH

      CCTK_INT,INTENT(IN) :: 
     &     conformal_state
      INTEGER,INTENT(IN) :: 
     &     igrid,l,m,Psi_power,myproc
      CCTK_INT,INTENT(IN) ::
     &     Nt,Np,nx,ny,nz,all_modes,do_momentum,do_spin,interpolation_order
      INTEGER,INTENT(IN) ::
     &     do_ADMmass(2)
      CCTK_REAL,INTENT(IN) :: 
     &     origin(3),Dx,Dy,Dz,eta
      CCTK_REAL,INTENT(IN),DIMENSION(nx,ny,nz) :: 
     &     Psi,g00,gxx,gxy,gxz,gyy,gyz,gzz,
     &     hxx,hxy,hxz,hyy,hyz,hzz
      CCTK_REAL,INTENT(IN),DIMENSION(nx) :: x
      CCTK_REAL,INTENT(IN),DIMENSION(ny) :: y
      CCTK_REAL,INTENT(IN),DIMENSION(nz) :: z
      CCTK_REAL,INTENT(INOUT),DIMENSION(nx,ny,nz) ::
     &     Extract_Temp3d
      CCTK_STRING,INTENT(IN) ::
     &     interpolation_operator

      CCTK_REAL,INTENT(OUT) :: 
     &     ADMmass(2),mass,rsch,Qodd(:,2:,0:),Qeven(:,2:,0:),dtaudt,
     &     momentum(3),spin(3)


      END SUBROUTINE
         
      END INTERFACE

      END MODULE D3_extract_int
