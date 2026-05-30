
#include "cctk.h"

      MODULE D2_extract_int

c     ------------------------------------------------------------------

      INTERFACE 

      SUBROUTINE D2_extract(eta,igrid,Nt,Np,theta,phi,all_modes,l,m,
     &        g00,g11,g12,g13,g22,g23,g33,dg22,dg23,dg33,
     &        do_ADMmass,ADMmass_int1,ADMmass_int2,
     &        do_momentum,momentum_int1,momentum_int2,momentum_int3,
     &        do_spin,spin_int1,spin_int2,spin_int3,
     &        mass,rsch,Qodd,Qeven,ADMmass,momentum,spin,dtaudt)

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: 
     &     igrid,l,m
      CCTK_INT,INTENT(IN) ::
     &     Nt,Np,all_modes,do_momentum,do_spin
      INTEGER,INTENT(IN) :: do_ADMmass(2)
      CCTK_REAL,INTENT(IN) :: eta
      CCTK_REAL,INTENT(IN),DIMENSION (Nt) :: theta
      CCTK_REAL,INTENT(IN),DIMENSION (Np) :: phi
      CCTK_REAL,INTENT(IN),DIMENSION (Nt,Np) :: 
     &     g00,g12,g13,g23,dg23,
     &     ADMmass_int1, ADMmass_int2,
     &     momentum_int1,momentum_int2,momentum_int3,
     &     spin_int1,spin_int2,spin_int3

      CCTK_REAL,INTENT(INOUT),DIMENSION (Nt,Np) ::
     &     g11,g22,g33,dg22,dg33

      CCTK_REAL,INTENT(OUT) :: 
     &     ADMmass(2),mass,rsch,Qodd(:,2:,0:),Qeven(:,2:,0:),dtaudt,
     &     momentum(3),spin(3)
      END SUBROUTINE
      
      END INTERFACE
      
      END MODULE D2_extract_int


