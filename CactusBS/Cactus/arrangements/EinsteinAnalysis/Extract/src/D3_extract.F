
#include "cctk.h"

c     ==================================================================

      SUBROUTINE D3_extract(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,igrid,
     &     origin,myproc,interpolation_operator,interpolation_order,Nt,Np,
     &     nx,ny,nz,all_modes,
     &     l,m,x,y,z,Dx,Dy,Dz,Psi_power,Psi,
     &     g00,gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &     eta,ADMmass,momentum,spin,mass,rsch,Qodd,Qeven,Extract_temp3d,dtaudt)

c     ------------------------------------------------------------------
c
c     Extract odd and even parity gauge invariant variables on a 
c     2-surface defined by constant isotropic radius ETA from ORIGIN,
c     also returns the corresponding Schwarzschild radius RSCH.
c   
c     Works for a full grid (IGRID=0) or an octant (IGRID=1), given 
c     non-conformal (PSI=1,Gij) or conformal (PSI^PSI_POWER,Gij).
c
c     IMPORTANT : Nt and Np must both be even numbers (so that Simpsons
c                 rule works later)
c
c     Variables in:
c     ------------
c     igrid     Grid type (0=full,1=octant)
c               [Could extract in an octant from a full grid, but why
c                would you want to]
c     origin    The position of the origin of extraction in the x,y,z
c               coordinate system
c     myproc    Cactus variable, must be zero if not using Cactus
c     Nt,Np     Number of theta,phi grid DIVISIONS to use [even numbers]
c     all_modes if true the extract all l,m modes up to l
c     l,m       Spherical harmonic to extract if ALL_MODES is false. 
c               otherwise m is ignored and l is the maximum l mode
c     x,y,z     Cartesian coordinates in the 3-space, hopefully 
c               isotropic about ORIGIN
c     Dx,Dy,Dz  Grid spacings
c     Psi_power The power of Psi which is being passed in, usually one
c               but may be four
c     Psi       The conformal factor, set it to one if non-conformal 
c               data is passed in
c     gij       The 3-metric components in cartesian coordinates
c     eta       The coordinate radius from ORIGIN for extraction, must
c               obviously give a 2-surface which is contained in the
c               grid
c
c     Variables out:
c     -------------
c     ADMmass   Estimate of ADM mass 
c     momentum  Estimate of momentum
c     spin      Estimate of spin
c     mass      The extracted mass
c     rsch      The extracted Schwarzschild radius
c     Qodd      The extracted odd parity gauge invariant variable
c     Qeven     The extracted even parity gauge invariant variable     
c     
c     ------------------------------------------------------------------

      USE D3_to_D2_int
      USE cartesian_to_spherical_int
      USE unphysical_to_physical_int
      USE D2_extract_int

c     ------------------------------------------------------------------

      IMPLICIT NONE

c     Input variables

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
     &     Extract_temp3d
      CCTK_STRING,INTENT(IN) ::
     &     interpolation_operator

c     Output variables
      CCTK_REAL,INTENT(OUT) :: 
     &     ADMmass(2),mass,rsch,Qodd(:,2:,0:),Qeven(:,2:,0:),dtaudt,
     &     momentum(3),spin(3)

c     Local variables passed on
      CCTK_REAL :: 
     &    theta(Nt),phi(Np)
      CCTK_REAL,DIMENSION(Nt,Np) :: 
     &    Psis,g00s,gxxs,gxys,gxzs,gyys,gyzs,gzzs,dPsis,dgxxs,dgxys,dgxzs,
     &    dgyys,dgyzs,dgzzs,grr,grt,grp,gtt,gtp,gpp,dgtt,dgtp,dgpp,
     &    ADMmass_int1,ADMmass_int2,
     &    momentum_int1,momentum_int2,momentum_int3,
     &    spin_int1,spin_int2,spin_int3

c     Local variables here only
      INTEGER :: 
     &    i,j
      CCTK_REAL,PARAMETER ::
     &    half = 0.5D0,
     &    one  = 1.0D0,
     &    two  = 2.0D0
      CCTK_REAL :: 
     &    Pi,Dt,Dp

c     ------------------------------------------------------------------


c     ------------------------------------------------------------------
c
c     1. Specify polar coordinates on chosen 2-surface
c
c     ------------------------------------------------------------------

c      print*,'Entering D3_extract'
c      print*,nx,ny,nz
c      print*,'D3_extract Qodd'
c      print*,lbound(Qodd,1),lbound(Qodd,2),lbound(Qodd,3)
c      print*,size(Qodd,1),size(Qodd,2),size(Qodd,3)
      Pi = two*ASIN(one)

c     Set polar gridspacing

      SELECT CASE (igrid)
      CASE (0)                  ! full grid
        Dt = Pi/DBLE(Nt)
        Dp = two*Pi/DBLE(Np)
      CASE (1)                  ! octant grid
        Dt = half*Pi/DBLE(Nt)
        Dp = half*Pi/DBLE(Np)
      CASE (2)                  ! cartoon
        Dt = Pi/DBLE(Nt)
        Dp = 2D0*Pi
      CASE (3)                  ! bitant
        Dt = half*Pi/DBLE(Nt)
        Dp = two*Pi/DBLE(Np)
      CASE DEFAULT
        WRITE (*,*) "Grid incorrectly set in D3_extract"
        STOP
      END SELECT
         

c     Set coordinates 

      DO i = 1, Nt
         theta(i) = (DBLE(i) - half)* Dt
      ENDDO
      DO j = 1, Np
         phi(j) = (DBLE(j) - half)* Dp
      ENDDO

      IF (igrid == 2) phi(1) = 0D0
      

c     ------------------------------------------------------------------
c
c     2. Project quantities onto the 2-surface
c
c     ------------------------------------------------------------------

c      print*,nx,ny,nz,Nt,Np
c      print*,'Calling D3_to_D2'
      CALL D3_to_D2(cctkGH,conformal_state,do_ADMmass,do_momentum,do_spin,
     &     Psi_power,origin,myproc,interpolation_operator,interpolation_order,
     &     Dx,Dy,Dz,Psi,
     &     g00,gxx,gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &     x,y,z,eta,Nt,Np,nx,ny,nz,theta,phi,Psis,g00s,gxxs,gxys,gxzs,
     &     gyys,gyzs,gzzs,dPsis,dgxxs,dgxys,dgxzs,dgyys,dgyzs,dgzzs,
     &     ADMmass_int1,ADMmass_int2,momentum_int1,momentum_int2,
     &     momentum_int3,spin_int1,spin_int2,spin_int3,Extract_temp3d)
c      print*,'Exited D3_to_D2'


      if (myproc .eq. 0) then 

c     ------------------------------------------------------------------
c
c     3. Transform tensor components to polar coordinates on 2-surface
c
c     ------------------------------------------------------------------

         CALL cartesian_to_spherical(theta,phi,eta,gxxs,gxys,gxzs,gyys,
     &        gyzs,gzzs,grr,grt,grp,gtt,gtp,gpp,dgxxs,dgxys,dgxzs,dgyys,
     &        dgyzs,dgzzs,dgtt,dgtp,dgpp)


c     ------------------------------------------------------------------
c
c     4. Calculate physical quantities on 2-surface
c
c     ------------------------------------------------------------------

         if (conformal_state > 0) then
           CALL unphysical_to_physical(grr,grt,grp,gtt,gtp,gpp,dgtt,dgtp,
     &        dgpp,Psis,dPsis,Psi_power)
         end if

c     ------------------------------------------------------------------
c
c     5. Extract gauge invariant quantities on 2-surface
c
c     ------------------------------------------------------------------

c      print*,'Calling D2_extract'
         CALL D2_extract(eta,igrid,Nt,Np,theta,phi,all_modes,l,m,g00s,grr,
     &        grt,grp,gtt,gtp,gpp,dgtt,dgtp,dgpp,
     &        do_ADMmass,ADMmass_int1,ADMmass_int2,
     &        do_momentum,momentum_int1,momentum_int2,momentum_int3,
     &        do_spin,spin_int1,spin_int2,spin_int3,
     &        mass,rsch,Qodd,Qeven,ADMmass,momentum,spin,dtaudt)
c      print*,Qodd
c      print*,Qeven
c      print*,'Exiting D2_extract'
      

c     End of myproc = 0 conditional
      endif

      END SUBROUTINE D3_extract



