 /*@@
   @file      WaveToy.F
   @date      
   @author    Tom Goodale
   @desc 
              Evolution routines for the wave equation solver
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    WaveToyF90_Evolution
   @date       
   @author     Tom Goodale
   @desc 
               Evolution for the wave equation
   @enddesc 
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyf90_Evolution(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

!     Declare local variables
      INTEGER  i,j,k
      INTEGER  istart, jstart, kstart, iend, jend, kend

      CCTK_REAL  dx,dy,dz,dt
      CCTK_REAL  dx2,dy2,dz2,dt2
      CCTK_REAL  dx2i,dy2i,dz2i

      CCTK_REAL  factor

	  !------------------------
	  
	  CCTK_REAL M0,L0,K1,L1,K2a,L2a,K2b,L2b,K2c,L2c
	  CCTK_REAL K2d,L2d,K3a,L3a,K3b,L3b,K3c,L3c,K3d,L3d,K4a
	  CCTK_REAL L4a,K4b,L4b,K4c,L4c,K4d,L4d,K4e,L4e,K4f,L4f
	  CCTK_REAL K4g,L4g,K4h,L4h,K4i,L4i
	  CCTK_REAL M,M1,aa1,qq,a1,a2,m_1,m_2,aa
          CCTK_REAL Z1,Z2,Spar,risco,Lisco,Eisco,aeff,dm,Dpar,eta,one
          one=1.0d0
	  
          print*,x(0,0,0),y(0,0,0),z(0,0,0)	  
	  
c Some known constants 
	  M0  =  0.951507
	  L0  =  0.686710
	  K1  = -0.051379
	  L1  =  0.613247
	  K2a = -0.004804 
	  L2a = -0.145427
	  K2b = -0.054522
	  L2b = -0.115689
	  K2c = -0.000022
	  L2c = -0.005254
	  K2d =  1.995246
	  L2d =  0.801838
	  K3a =  0.007064
	  L3a = -0.073839
	  K3b = -0.017599
	  L3b =  0.004759
	  K3c = -0.119175 
	  L3c = -0.078377
	  K3d =  0.025000
	  L3d =  1.585809
	  K4a = -0.068981 
	  L4a = -0.003050
	  K4b = -0.011383
	  L4b = -0.002968
	  K4c = -0.002284
	  L4c =  0.004364
	  K4d = -0.165658
	  L4d = -0.047204
	  K4e =  0.019403
	  L4e = -0.053099 
	  K4f =  2.980990
	  L4f =  0.953458
	  K4g =  0.020250
	  L4g = -0.067998
	  K4h = -0.004091 
	  L4h =  0.001629
	  K4i =  0.078441
	  L4i = -0.066693

      M=M_total
      M1=Final_Mass
      aa1=Final_Spin
	  
       print *, 'are we here?'	  

!     Set up shorthands
!     -----------------
      dx = CCTK_DELTA_SPACE(1)
      dy = CCTK_DELTA_SPACE(2)
      dz = CCTK_DELTA_SPACE(3)
      dt = CCTK_DELTA_TIME

      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz
      dt2 = dt*dt

      dx2i = 1.0/dx2
      dy2i = 1.0/dy2
      dz2i = 1.0/dz2

      istart = 2
      jstart = 2
      kstart = 2
      
      iend = cctk_lsh(1)-1
      jend = cctk_lsh(2)-1
      kend = cctk_lsh(3)-1
      print*,iend
      factor = 2*(1 - (dt2)*(dx2i + dy2i + dz2i))
      open(1, file = 'brute.dat', status = 'new')  

!     Do the evolution
!     ----------------
      do k = kstart, kend
         do j = jstart, jend
            do i = istart, iend
		
                print*, 'IN THE LOOP'		
		a1=y(i,j,k)
		a2=z(i,j,k)
		m_1 = M * x(i,j,k) / ( 1.0 + x(i,j,k) ) 
		m_2 = M / ( 1.0 + x(i,j,k) )
                print*,a1,a2,x(i,j,k)	
		dm = ( m_1 - m_2 ) / M
		eta = m_1 * m_2 / M**2
		Spar = ( a1 * m_1**2 + a2 * m_2**2 ) / m**2
                aeff=Spar
                Dpar = (a2*m_2-a1*m_1)/m
                Z1 = 1+(1-aeff**2)**(1/3)*((1+aeff)**(1/3)+(1-aeff)**(1/3))
                Z2 = sqrt(3*aeff**2+Z1**2)
                risco = 3+Z2-SIGN(one,aeff)*sqrt((3-Z1)*(3+Z1+2*Z2))
                Lisco = 2/sqrt(3*risco)*(3*sqrt(risco)-2*aeff)
                Eisco = (1-2/risco+aeff/risco**(1.5))/sqrt(1-3/risco+2*aeff/risco**(1.5))
				
                M1=(4*eta)**2*(M0 + K1*Spar + K2a*Dpar*dm + K2b*Spar**2 + K2c*Dpar**2 + 
     1		K2d*dm**2 + K3a*Dpar*Spar*dm + K3b*Spar*Dpar**2 + K3c*Spar**3 + K3d*Spar*dm**2 +
     2          K4a*Dpar*Spar**2*dm + K4b*Dpar**3*dm + K4c*Dpar**4 + K4d*Spar**4 +
     3          K4e*Dpar**2*Spar**2 + K4f*dm**4 + K4g*Dpar*dm**3 + K4h*Dpar**2*dm**2 + K4i*Spar**2*dm**2) + (1+eta*(Eisco+11))*dm**6
    
	        aa=(4*eta)**2*(L0 + L1*Spar+L2a*Dpar*dm+L2b*Spar**2+L2c*Dpar**2+L2d*dm**2 +
     1		L3a*Dpar*Spar*dm+L3b*Spar*Dpar**2+L3c*Spar**3+L3d*Spar*dm**2+L4a*Dpar*Spar**2*dm+
     2          L4b*Dpar**3*dm+L4c*Dpar**4+L4d*Spar**4+L4e*Dpar**2*Spar**2+L4f*dm**4+
     3          L4g*Dpar*dm**3+L4h*Dpar**2*dm**2+L4i*Spar**2*dm**2) + Spar*(1+8*eta)*dm**4+eta*Lisco*dm**6
                IF ( abs((aa-aa1)/aa1) .LE. precision .AND. abs((M-M1)/M1) .LE. precision ) THEN 
                                    print*,'saving'
				       write(1,*) x(i,j,k),y(i,j,k),z(i,j,k)   

				 
                END IF 
          
               
               phi(i,j,k) = factor*phi_p(i,j,k) -
     1              phi_p_p(i,j,k) + (dt2) * 
     2              ((phi_p(i+1,j,k)+phi_p(i-1,j,k))*dx2i
     3              +(phi_p(i,j+1,k)+phi_p(i,j-1,k))*dy2i
     4              +(phi_p(i,j,k+1)+phi_p(i,j,k-1))*dz2i)

            end do
         end do
      end do

                close(1)         
      end subroutine WaveToyF90_Evolution


 /*@@
   @routine    WaveToyF90_Boundaries
   @date       
   @author     Tom Goodale
   @desc 
               Boundary conditions for the wave equation
   @enddesc 
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyF90_Boundaries(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_INT, parameter :: izero = 0
      integer, parameter :: ik = kind(izero)

c     Local declarations
      CCTK_INT :: ierr

      CHARACTER (len=100) :: boundary
      INTEGER :: length


      ierr = 0

c     The "bound" parameter needs to be converted into a Fortran string.
      call CCTK_FortranString(length,bound,boundary)

c     Apply the outer boundary conditions
c     -----------------------------------
c     Note: In each of the following calls to Boundary_SelectVarForBC,
c     default arguments are used, so an invalid table handle of -1 can
c     be passed
      if (CCTK_EQUALS(bound,"flat") .or. CCTK_EQUALS(bound,"static") .or.
     $    CCTK_EQUALS(bound,"radiation") .or. CCTK_EQUALS(bound,"robin") .or.
     $    CCTK_EQUALS(bound,"none") ) then
         ierr = Boundary_SelectVarForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "wavetoy::phi", boundary)
      else if (CCTK_EQUALS(bound,"zero")) then
         ierr = Boundary_SelectVarForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "wavetoy::phi", "Scalar")
      end if

      if (ierr < 0)  then
        call CCTK_WARN(0,"WaveToyF90_Boundaries: Error in selecting boundary condition")
      end if

      end subroutine WaveToyF90_Boundaries
