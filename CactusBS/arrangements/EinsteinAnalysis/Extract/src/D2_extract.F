
#include "cctk.h"

c     ========================================================================

      SUBROUTINE D2_extract(eta,igrid,Nt,Np,theta,phi,all_modes,l,m,
     &     g00,g11,g12,g13,g22,g23,g33,dg22,dg23,dg33,
     &     do_ADMmass,ADMmass_int1,ADMmass_int2,
     &     do_momentum,momentum_int1,momentum_int2,momentum_int3,
     &     do_spin,spin_int1,spin_int2,spin_int3,
     &     mass,rsch,Qodd,Qeven,ADMmass,momentum,spin,dtaudt)

c     ------------------------------------------------------------------------
c
c     Extract the odd and even parity gauge invariant variables
c     at a given isotropic radius (eta) from the center of a 
c     slightly perturbed Schwarzschild spacetime.
c     
c     The polar coordinates and metric variables must be given on 
c     an equally spaced grid, which is offset from the axis.
c
c     Notice that for octant symmetry we only need calculate the real 
c     parts of the gauge invariant variables for the even-l,m modes.
c     This subroutine is optimised (!!) for octant symmetry.
c
c     ------------------------------------------------------------------------

      IMPLICIT NONE

c     Input variables
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
     &     ADMmass_int1,ADMmass_int2,
     &     momentum_int1,momentum_int2,momentum_int3,
     &     spin_int1,spin_int2,spin_int3

      CCTK_REAL,INTENT(INOUT),DIMENSION (Nt,Np) ::
     &     g11,g22,g33,dg22,dg33

c     Output variables
      CCTK_REAL,INTENT(OUT) :: 
     &     ADMmass(2),mass,rsch,Qodd(:,2:,0:),Qeven(:,2:,0:),dtaudt,
     &     momentum(3),spin(3)

c     Local Variables
      LOGICAL, PARAMETER :: 
     &    debug = .FALSE.
      INTEGER ::  
     &    i,j,il,im,lmin,lmax,mmin,mmax,lstep,mstep
      CCTK_REAL,PARAMETER ::
     &    zero  = 0.0D0,
     &    half  = 0.5D0,
     &    one   = 1.0D0,
     &    two   = 2.0D0,
     &    four  = 4.0D0,
     &    six   = 6.0D0,
     &    eight = 8.0D0
      CCTK_REAL :: 
     &     Dt,Dp,st,ist,drschdri,dridrsch,
     &     g22comb,g23comb,g33comb,S,rsch2,Lambda,
     &     Pi,
     &     fac_h1,fac_H2,fac_K,fac_G,fac_c1,fac_c2,
     &     fac_dG,fac_dK,fac_dc2,
     &     sph_g11,sph_g22,sph_g33,sph_dg22,
     &     sphere_int
      CCTK_REAL,DIMENSION(2) :: 
     &     Y,Y1,Y2,Y3,Y4,h1,H2,K,G,c1,c2,dG,dK,dc2
      CCTK_REAL,DIMENSION(Nt+1,Np+1) :: 
     &     temp,h1i,H2i,Gi,Ki,c1i,c2i,dGi,dKi,dc2i


c     ------------------------------------------------------------------------
c
c     0. Initial things
c
c     ------------------------------------------------------------------------

c      print*,'D2_extract Qodd'
c      print*,lbound(Qodd,1),lbound(Qodd,2),lbound(Qodd,3)
c      print*,size(Qodd,1),size(Qodd,2),size(Qodd,3)
      Pi    = two*ASIN(one)

      SELECT CASE (igrid)
      CASE (0)                  ! full grid
        Dt = Pi/DBLE(Nt)
        Dp = two*Pi/DBLE(Np)
        lstep = 1
        mstep = 1
      CASE (1)                  ! octant grid
        Dt = half*Pi/DBLE(Nt)
        Dp = half*Pi/DBLE(Np)
        lstep = 2
        mstep = 2
      CASE (2)                  ! cartoon
        Dt = Pi/DBLE(Nt)
        Dp = 2D0*Pi
        lstep = 2
        mstep = 2
      CASE (3)                  ! bitant
        Dt = half*Pi/DBLE(Nt)
        Dp = two*Pi/DBLE(Np)         
        lstep = 1
        mstep = 1
      CASE DEFAULT
        WRITE (*,*) "Grid incorrectly set in D2_extract"
        STOP
      END SELECT
         
c     Calculate ADM mass
c     ------------------
      IF (do_ADMmass(1) == 1) THEN

c      print*,'1'
c       Calculate ADM mass using standard formula 
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = ADMmass_int1(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         ADMmass(1) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)

      ENDIF
      IF (do_ADMmass(2) == 1) THEN
 
c      print*,'2'
c       Calculate ADM mass using conformal formula
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = ADMmass_int2(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         ADMmass(2) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
         
      END IF

c     Calculate momentum
c     ------------------
      IF (do_momentum == 1) THEN

c      print*,'3'
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = momentum_int1(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         momentum(1) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
         
c      print*,'4'
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = momentum_int2(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         momentum(2) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
 
c      print*,'5'
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = momentum_int3(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         momentum(3) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
 
      ENDIF

c     Calculate spin
c     --------------
      IF (do_spin == 1) THEN

c      print*,'6'
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = spin_int1(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         spin(1) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
         
c      print*,'7'
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = spin_int2(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         spin(2) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
 
c      print*,'8'
         DO j = 2, Np+1         
            DO i = 2, Nt+1
               temp(i,j) = spin_int3(i-1,j-1)*SIN(theta(i-1))
            ENDDO
         ENDDO
         spin(3) = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
 
      ENDIF

c     ------------------------------------------------------------------
c
c     1. Find Spherical Parts of Metric Components
c
c     ------------------------------------------------------------------

c     ... g00
c      print*,'9'
      DO j = 2, Np+1
        DO i = 2, Nt+1
          temp(i,j) = g00(i-1,j-1)*SIN(theta(i-1))
        ENDDO
      ENDDO
      dtaudt = sqrt(one/(four*Pi)*sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp))

c     ... g11
c      print*,'10'
      DO j = 2, Np+1
         DO i = 2, Nt+1
            temp(i,j) = g11(i-1,j-1)*SIN(theta(i-1))
         ENDDO
      ENDDO
      sph_g11 = one/(four*Pi)*sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)

c     ... g22
c      print*,'11'
      DO j = 2, Np+1
         DO i = 2, Nt+1
            temp(i,j) = g22(i-1,j-1)*SIN(theta(i-1))
         ENDDO
      ENDDO
      sph_g22 = one/(four*Pi)*sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)

c     ... dg22
c      print*,'12'
      DO j = 2, Np+1
         DO i = 2, Nt+1
            temp(i,j) = dg22(i-1,j-1)*SIN(theta(i-1))
         ENDDO
      ENDDO
      sph_dg22 = one/(four*Pi)*sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)


c     ... g33/sin(theta)**2
c      print*,'13'
      DO j = 2, Np+1
         DO i = 2, Nt+1
            temp(i,j) = g33(i-1,j-1)/SIN(theta(i-1))
         ENDDO
      ENDDO
      sph_g33 = one/(four*Pi)*sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)


c     Calculate error in orthonormality of spherical harmonic
c      DO j = 2, Np+1
c         DO i = 2, Nt+1
c            CALL spher_harm_combs(theta(i),phi(j),l,m,Y,Y1,Y2,Y3,Y4)
c            temp(i,j) = (Y(1)*Y(1)+Y(2)*Y(2))*sin(theta(i-1))
c         ENDDO
c      ENDDO
c      error = sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c      print *,"Error is",one-error



c     ------------------------------------------------------------------
c
c     1. Compute the Schwarzschild radius
c   
c     ------------------------------------------------------------------

c      print*,'14'
      DO j = 2, Np+1
         DO i = 2, Nt+1
c            temp(i,j) = g22(i-1,j-1)*SIN(theta(i-1))
            temp(i,j) = (g22(i-1,j-1)+g33(i-1,j-1)/SIN(theta(i-1))**2)
     &           *SIN(theta(i-1))/(two)
         ENDDO
      ENDDO

      rsch2= one/(four*Pi)*sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)

      rsch = SQRT(rsch2)


c     ------------------------------------------------------------------
c
c     2. Compute the derivative of the schwarzschild radius with  
c        respect to the isotropic radius eta.
c
c     ------------------------------------------------------------------

c      print*,'15'
      DO j = 2, Np+1         
         DO i = 2, Nt+1
c           temp(i,j) = (dg22(i-1,j-1))*SIN(theta(i-1))
           temp(i,j) = (dg22(i-1,j-1)+dg33(i-1,j-1)/SIN(theta(i-1))**2)
     &             *SIN(theta(i-1))/(two)
         ENDDO
      ENDDO

      drschdri = one/(eight*Pi*rsch)*
     &     sphere_int(temp,igrid,int(Nt+1),int(Np+1),Dt,Dp)

      dridrsch = one/drschdri


c     ------------------------------------------------------------------
c
c     3. Calculate the Schwarzschild mass and S factor
c
c     ------------------------------------------------------------------

      S = drschdri**2/sph_g11
      mass = rsch*(one-S)/two



c     ------------------------------------------------------------------
c
c     4. Subtract spherical part from metric
c        Do not know how important this is
c
c     ------------------------------------------------------------------

c      print*,'16'
      DO i=1,Nt
         DO j=1,Np
            st  = SIN(theta(i))
            g11(i,j) = g11(i,j) - sph_g11
            g22(i,j) = g22(i,j) - sph_g22
            g33(i,j) = g33(i,j) - sph_g22*st**2
            dg22(i,j) = dg22(i,j) - sph_dg22
            dg33(i,j) = dg33(i,j) - sph_dg22*st**2
         ENDDO
      ENDDO


      IF (all_modes == 0) THEN
         lmin = l ; lmax = l
      ELSE
         lmin = 2 ; lmax = l
      ENDIF
      
      loop_l: DO il = lmin,lmax,lstep
         
         IF (all_modes == 0) THEN
            mmin = m ; mmax = m
         ELSE
            mmin = 0 ; mmax = il
         END IF
         
         loop_m: DO im = mmin,mmax,mstep

c     ------------------------------------------------------------------
c     
c     5. Calculate all Regge-Wheeler variables
c
c     ------------------------------------------------------------------

c           Factors independent of angular coordinates
c           ------------------------------------------
            fac_h1  = dridrsch/DBLE(il*(il+1))
            fac_H2  = S*dridrsch**2
            fac_G   = one/(rsch**2*DBLE(il*(il+1)*(il-1)*(il+2)))
            fac_K   = half/rsch**2
            fac_c1  = dridrsch/DBLE(il*(il+1))
            fac_c2  = two/DBLE(il*(il+1)*(il-1)*(il+2))
            fac_dG  = fac_G
            fac_dK  = fac_K
            fac_dc2 = fac_c2
  
          
c     ------------------------------------------------------------------
c   
c     5a. Real parts of the Regge-Wheeler variables
c
c     ------------------------------------------------------------------

c           Calculate integrands

c      print*,'17'
            loop_phi1: DO j = 1, Np
               loop_theta1: DO i = 1, Nt

                  st  = SIN(theta(i))
                  ist = one/st
                  
c                  print*,'Calling spher_harm_coms',i,j
                  CALL spher_harm_combs(theta(i),phi(j),il,im,Y,Y1,Y2,Y3
     &                 ,Y4)
                  
                  h1i(i+1,j+1) =  st*g12(i,j)*Y1(1)+ist*g13(i,j)*Y2(1)    
                  H2i(i+1,j+1) =  st*g11(i,j)*Y(1)
                  Gi(i+1,j+1)  =  (st*g22(i,j)-ist*g33(i,j))*Y3(1)
     &                 +four*ist*g23(i,j)*Y4(1)
                  Ki(i+1,j+1) =  (st*g22(i,j)+ist*g33(i,j))*Y(1)              
                  c1i(i+1,j+1) =  g13(i,j)*Y1(1)-g12(i,j)*Y2(1)           
                  c2i(i+1,j+1) =  (g22(i,j)-ist**2*g33(i,j))*Y4(1)
     &                 -g23(i,j)*Y3(1)
                  
                  g22comb = dridrsch*dg22(i,j)-two/rsch*g22(i,j)
                  g23comb = dridrsch*dg23(i,j)-two/rsch*g23(i,j)
                  g33comb = dridrsch*dg33(i,j)-two/rsch*g33(i,j)
                  
                  dGi(i+1,j+1) =  (st*g22comb-ist*g33comb)*Y3(1)
     &                 +four*ist*g23comb*Y4(1)                  
                  dKi(i+1,j+1) =  (st*g22comb+ist*g33comb)*Y(1)
                  dc2i(i+1,j+1) =  (dg22(i,j)-ist**2*dg33(i,j))*Y4(1)
     &                 -dg23(i,j)*Y3(1)

               END DO loop_theta1    
            END DO loop_phi1
c            print*,'End loop'

c           Integrations over the 2-sphere

            h1(1) = fac_h1*sphere_int(h1i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'h1 done'
            H2(1) = fac_H2*sphere_int(H2i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'H2 done'
            G(1) = fac_G*sphere_int(Gi,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'G done'
            K(1) = fac_K*sphere_int(Ki,igrid,int(Nt+1),int(Np+1),Dt,Dp)
     &           +DBLE(il*(il+1))*half*G(1)
c            print*,'K done'
            c1(1) = fac_c1*sphere_int(c1i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'c1 done'
            c2(1) = fac_c2*sphere_int(c2i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'c2 done'
            dG(1) = fac_dG*sphere_int(dGi,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'dG done'
            dK(1) = fac_dK*sphere_int(dKi,igrid,int(Nt+1),int(Np+1),Dt,Dp)
     &           +DBLE(il*(il+1))*half*dG(1)
c            print*,'dK done'
            dc2(1) = fac_dc2*sphere_int(dc2i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
c            print*,'dc2 done'


c     ------------------------------------------------------------------
c   
c     5b. Imaginary parts of the Regge-Wheeler variables
c
c     ------------------------------------------------------------------

c      print*,'17 1/2'
            complex_parts: IF (igrid /= 1 .AND. igrid /=2) THEN

c              Calculate integrands

c      print*,'18'
               loop_phi2: DO j = 1, Np

                  loop_theta2: DO i = 1, Nt
                  
                     st  = SIN(theta(i))
                     ist = one/st
                  
                     CALL spher_harm_combs(theta(i),phi(j),il,im,Y,Y1,Y2
     &                    ,Y3,Y4)
                  
                     h1i(i+1,j+1) =-st*g12(i,j)*Y1(2)-ist*g13(i,j)*Y2(2)
                     H2i(i+1,j+1) = -st*g11(i,j)*Y(2) 
                     Gi(i+1,j+1)  = -(st*g22(i,j)-ist*g33(i,j))*Y3(2)
     &                    -four*ist*g23(i,j)*Y4(2)                  
                     Ki(i+1,j+1) = -(st*g22(i,j)+ist*g33(i,j))*Y(2)   
                     c1i(i+1,j+1) = -g13(i,j)*Y1(2)+g12(i,j)*Y2(2) 
                     c2i(i+1,j+1) = -(g22(i,j)-ist**2*g33(i,j))*Y4(2)
     &                    +g23(i,j)*Y3(2)
                     
                     g22comb = dridrsch*dg22(i,j)-two/rsch*g22(i,j)
                     g23comb = dridrsch*dg23(i,j)-two/rsch*g23(i,j)
                     g33comb = dridrsch*dg33(i,j)-two/rsch*g33(i,j)
                     
                     dGi(i+1,j+1) = -(st*g22comb-ist*g33comb)*Y3(2)
     &                    -four*ist*g23comb*Y4(2) 
                     dKi(i+1,j+1) = -(st*g22comb+ist*g33comb)*Y(2)  
                     dc2i(i+1,j+1) = -(dg22(i,j)-ist**2*dg33(i,j))*Y4(2)
     &                    +dg23(i,j)*Y3(2)
                  
                  END DO loop_theta2    
               END DO loop_phi2


c              Integrate over 2-sphere

               h1(2) = fac_h1*sphere_int(h1i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
               H2(2) = fac_H2*sphere_int(H2i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
               G(2) = fac_G*sphere_int(Gi,igrid,int(Nt+1),int(Np+1),Dt,Dp)
               K(2) = fac_K*sphere_int(Ki,igrid,int(Nt+1),int(Np+1),Dt,Dp)
     &              +DBLE(il*(il+1))*half*G(2)
               c1(2) = fac_c1*sphere_int(c1i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
               c2(2) = fac_c2*sphere_int(c2i,igrid,int(Nt+1),int(Np+1),Dt,Dp)
               dG(2) = fac_dG*sphere_int(dGi,igrid,int(Nt+1),int(Np+1),Dt,Dp)
               dK(2) = fac_dK*sphere_int(dKi,igrid,int(Nt+1),int(Np+1),Dt,Dp)
     &              +DBLE(il*(il+1))*half*dG(2)
               dc2(2) = fac_dc2*sphere_int(dc2i,igrid,int(Nt+1),int(Np+1),Dt,Dp)

            ELSE

               h1(2) = zero
               H2(2) = zero
               G(2) = zero
               K(2) = zero
               c1(2) = zero
               c2(2) = zero
               dG(2) = zero
               dK(2) = zero
               dc2(2) = zero
               

            END IF complex_parts
            


c     ------------------------------------------------------------------
c     
c     6. Construct gauge invariant variables
c
c     ------------------------------------------------------------------

c            print*,'Qodd',il,im
c            print*,lbound(Qodd,1),lbound(Qodd,2),lbound(Qodd,3)
c            print*,size(Qodd,1),size(Qodd,2),size(Qodd,3)
            Qodd(:,il,im) = SQRT(two*DBLE((il+2)*(il+1)*il*(il-1)))*S/
     &           rsch*(c1+half*(dc2-two/rsch*c2))
            
c            print*,'Lambda'
            Lambda = DBLE((il-1)*(il+2))+3.0D0*(one-S)
            
c            print*,'Qeven'
            Qeven(:,il,im) = one/Lambda*SQRT(two*DBLE((il-1)*(il+2))/
     &           DBLE(il*(il+1)))*(DBLE(il*(il+1))*S*(rsch**2*dG-
     &           two*h1)+two*rsch*S*(H2-rsch*dK)
     &           +Lambda*rsch*K)
            

         END DO loop_m
         
      END DO loop_l



c     ------------------------------------------------------------------
c     
c     7. Write some diagnostics if required (for last l,m to be used)
c
c     ------------------------------------------------------------------

      IF (debug) THEN
         WRITE(101,*) eta,h1
         WRITE(102,*) eta,H2
         WRITE(103,*) eta,G
         WRITE(104,*) eta,K
         WRITE(105,*) eta,c1
         WRITE(106,*) eta,c2
         WRITE(107,*) eta,dG
         WRITE(108,*) eta,dK
         WRITE(109,*) eta,dc2
         WRITE(201,*) eta,mass
         WRITE(202,*) eta,rsch
         WRITE(203,*) eta,drschdri
         WRITE(301,*) eta,Qeven(:,lmax,mmax)
         WRITE(302,*) eta,Qodd(:,lmax,mmax)
      ENDIF

c      print*,'Exiting D2_extract'

      END SUBROUTINE D2_extract

c     ===============================================================   
c
c
c
c
c
c     ===============================================================
      
      SUBROUTINE simpsons(n,Dx,f,intf)
      
c     Simpsons rule to evaluate 1D integral equation
      
c     Variables In:
c     n     number of points [must be odd]
c     Dx    step size
c     f(n)  integrand at each abscissa
      
c     Variables Out:
c     intf  evaluated integral
      
c     Variables Local:
c     i     index
      
c     ---------------------------------------------------------------
      
      IMPLICIT NONE
      
c     Input variables
      INTEGER,INTENT(IN) :: 
     &     n        
      CCTK_REAL,INTENT(IN) :: 
     &     Dx,f(n)      
      
c     Output variables
      CCTK_REAL,INTENT(OUT) :: 
     &     intf    
      
c     Local variables
      INTEGER :: 
     &     i
      CCTK_REAL,PARAMETER ::
     &     two   = 2.0D0,
     &     three = 3.0D0,
     &     four  = 4.0D0
      
c     ---------------------------------------------------------------

      IF (MOD(n,2) == 0) THEN
        WRITE(*,*) "Call to -simpsons- with even number of points"
        STOP
      END IF
            
      intf = f(1) + four*f(n-1) + f(n)
      
      DO i = 2, n-3,2
         
         intf = intf + four*f(i) + two*f(i+1)
         
      ENDDO
      
      intf = intf*Dx/three
      
      END SUBROUTINE simpsons
      
c     ===============================================================
c
c
c
c
c
c     ===============================================================

      FUNCTION sphere_int(temp,igrid,Nt,Np,Dt,Dp)

c     ---------------------------------------------------------------
c
c     Integration the function held in temp over the sphere. 
c     If the full grid is being used then Simpsons rule is used, 
c     if an octant is used then the symmetry means Simpsons rule 
c     becomes a simple sum.
c
c     ---------------------------------------------------------------

      IMPLICIT NONE

c     Input variables
      INTEGER,INTENT(IN) :: 
     &   Nt,Np,igrid
      CCTK_REAL,INTENT(IN) :: 
     &   temp(Nt,Np),Dt,Dp

c     Output variables
      CCTK_REAL ::
     &   sphere_int

c     Local variables
      INTEGER ::
     &   i,j
      CCTK_REAL :: 
     &   ft(Nt),fp(Np)
      CCTK_REAL,PARAMETER ::
     &   two  = 2.0D0,
     &   four = 4.0D0
   
c     ---------------------------------------------------------------  

      DO j = 2, Np

c        Integrate with respect to theta

         DO i = 2, Nt
            ft(i) = temp(i,j)
         ENDDO
         ft(1) = ft(2)     ! from symmetry across axis
         
         SELECT CASE (igrid)
         CASE (0)
            CALL simpsons(Nt,Dt,ft,fp(j))
         CASE (1)
            fp(j)=two*Dt*SUM(ft(2:Nt))
         CASE (2)
            CALL simpsons(Nt,Dt,ft,fp(j))
         CASE (3)
            fp(j)=two*Dt*SUM(ft(2:Nt))
         END SELECT

      ENDDO

c     Integrate with respect to phi

      fp(1) = fp(2)             ! from symmetry across axis


      SELECT CASE (igrid)
      CASE (0)
         CALL simpsons(Np,Dp,fp,sphere_int)
      CASE (1)
         sphere_int=four*Dp*SUM(fp(2:Np))
      CASE (2)
         sphere_int=Dp*fp(1)
      CASE (3)
         CALL simpsons(Np,Dp,fp,sphere_int)
      END SELECT   


      END FUNCTION sphere_int     

c     ==================================================================
c
c
c
c
c
c     ==================================================================

      SUBROUTINE spherical_harmonic(l,m,theta,phi,Ylm)

c     ------------------------------------------------------------------
c
c     Calculate the (l,m) spherical harmonic at given angular
c     coordinates. This number is in general complex, and
c
c      Ylm = Ylm(1) + i Ylm(2)
c
c      where
c     
c                a    ( 2 l + 1 (l-|m|)! )                      i m phi
c      Ylm = (-1) SQRT( ------- -------  ) P_l|m| (cos(theta)) e
c                     (   4 Pi  (l+|m|)! ) 
c
c      and where
c
c      a = m/2 (sign(m)+1)
c
c     ------------------------------------------------------------------

      IMPLICIT NONE

c     Input variables
      INTEGER ::
     &   l,m
      CCTK_REAL ::  
     &   theta,phi

c     Output variables
      CCTK_REAL :: 
     &   Ylm(2)

c     Local variables
      INTEGER ::  
     &    i
      CCTK_REAL,PARAMETER ::
     &    one  = 1.0D0,
     &    two  = 2.0D0,
     &    four = 4.0D0
      CCTK_REAL ::
     &    a,Pi,fac,plgndr

c     ------------------------------------------------------------------

      Pi = ACOS(-one)

      fac = one
      DO i = l-ABS(m)+1,l+ABS(m)
        fac = fac*DBLE(i)
      ENDDO
      fac = one/fac

c      a = (-one)**((m*ISIGN(m,1)/ABS(m)+m)/2)*SQRT(DBLE(2*l+1)
c     &    /four/Pi*fac)*plgndr(l,ABS(m),cos(theta))

      a = (-one)**MAX(m,0)*SQRT(DBLE(2*l+1)
     &    /four/Pi*fac)*plgndr(l,ABS(m),cos(theta))
      Ylm(1) = a*COS(DBLE(m)*phi)
      Ylm(2) = a*SIN(DBLE(m)*phi)

      END SUBROUTINE spherical_harmonic

c     ==================================================================
c
c
c
c
c
c     ==================================================================

      FUNCTION plgndr(l,m,x)

c     ------------------------------------------------------------------
c
c     From Numerical Recipes
c
c     Calculates the associated Legendre polynomial Plm(x).
c     Here m and l are integers satisfying 0 <= m <= l,
c     while x lies in the range -1 <= x <= 1
c
c     ------------------------------------------------------------------

c     Input variables
      INTEGER,INTENT(IN) :: 
     &    l,m
      CCTK_REAL,INTENT(IN) :: 
     &    x

c     Output variables
      CCTK_REAL :: 
     &    plgndr

c     Local Variables
      INTEGER :: 
     &    i,ll
      CCTK_REAL,PARAMETER ::
     &    one = 1.0D0,
     &    two = 2.0D0
      CCTK_REAL :: 
     &    pmm,somx2,fact,pmmp1,pll

c     ------------------------------------------------------------------

      pmm = one

      IF (m.GT.0) THEN
        somx2=SQRT((one-x)*(one+x))
        fact=one
        DO i=1,m
          pmm  = -pmm*fact*somx2
          fact = fact+two
        END DO
      ENDIF

      IF (l.EQ.m) THEN
         plgndr = pmm
      ELSE
         pmmp1 = x*(two*m+one)*pmm
         IF (l.EQ.m+1) THEN
            plgndr=pmmp1
         ELSE
            DO ll=m+2,l
               pll   = (x*DBLE(2*ll-1)*pmmp1-
     &              DBLE(ll+m-1)*pmm)/DBLE(ll-m)
               pmm   = pmmp1
               pmmp1 = pll
            END DO
            plgndr = pll
         ENDIF
      ENDIF

      END FUNCTION plgndr

c     ==================================================================
c
c
c
c
c
c     ==================================================================

      SUBROUTINE spher_harm_combs(theta,phi,l,m,Y,Y1,Y2,Y3,Y4)
c     ------------------------------------------------------------------
c
c     Calculates the various combinations of spherical harmonics needed
c     for the extraction (all are complex):     
c
c       Y  = Ylm
c       Y1 = Ylm,theta
c       Y2 = Ylm,phi
c       Y3 = Ylm,theta,theta-cot theta Ylm,theta-Ylm,phi,phi/sin^2 theta
c       Y4 = Ylm,theta,phi-cot theta Ylm,phi
c
c     The local variables Yplus is the spherical harmonic at (l+1,m)
c
c     ------------------------------------------------------------------

      IMPLICIT NONE

c     Input variables
      INTEGER ::
     &     l,m
      CCTK_REAL ::
     &     theta,phi

c     Output variables
      CCTK_REAL,DIMENSION(2) ::
     &     Y,Y1,Y2,Y3,Y4

c     Local variables
      INTEGER :: 
     &     i
      CCTK_REAL ::  
     &     Yplus(2),rl,rm,cot_theta,Pi
      CCTK_REAL,PARAMETER ::
     &     half = 0.5D0,
     &     one  = 1.0D0,
     &     two  = 2.0D0

c     ------------------------------------------------------------------

      Pi = ACOS(-one)

      rl = DBLE(l)
      rm = DBLE(m)

      cot_theta = COS(theta)/SIN(theta)

      CALL spherical_harmonic(l+1,m,theta,phi,Yplus)


c     Find Y ...

      CALL spherical_harmonic(l,m,theta,phi,Y)


c     Find Y1 ...

      DO i = 1,2
        Y1(i) = -(rl+one)*cot_theta*Y(i)+Yplus(i)/SIN(theta)*
     &          SQRT(((rl+one)**2-rm**2)*(rl+half)/(rl+one+half))
      ENDDO


c     Find Y2 ...

      Y2(1) = -rm*Y(2)
      Y2(2) =  rm*Y(1)


c     Find Y3 ...

      DO i = 1,2
        Y3(i) = -two*cot_theta*Y1(i)+(two*rm*rm/(SIN(theta)**2)
     &          -rl*(rl+one))*Y(i)
      ENDDO


c     Find Y4 ...

      Y4(1) = rm*(cot_theta*Y(2)-Y1(2))
      Y4(2) = rm*(Y1(1)-cot_theta*Y(1))


      END SUBROUTINE spher_harm_combs

c     ==================================================================




