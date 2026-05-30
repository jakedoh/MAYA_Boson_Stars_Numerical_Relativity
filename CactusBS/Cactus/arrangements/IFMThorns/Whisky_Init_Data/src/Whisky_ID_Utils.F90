 /*@@
   @file      ID_Util.F90
   @author    Tanja Bode
   @desc 
   Various utility functions.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"
#define PI  3.14159265358979323846  /* pi */

 /*@@
   @routine    RotateVector
   @date       Tue Jan 31 10:21:34 EST 2012
   @author     Tanja Bode 
   @desc 
   Helper routine to rotate initial v^i and B^i with respect to a shock
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RotateVector( vx, vy, vz, vparity, shock_direction, x, y, z )
  
  implicit none
  character(len=2), intent(in) :: shock_direction
  CCTK_INT, intent(in) :: vparity
  CCTK_REAL :: vx, vy, vz, x, y, z

  CCTK_REAL, dimension(3) :: ivec, rvec
  CCTK_REAL, dimension(3,3) :: Arot
  CCTK_REAL :: alpha, beta, gam !! Euler angles
  CCTK_INT :: i,j,euler_rotate

  !! Store initial vector
  ivec(1) = vx
  ivec(2) = vy
  ivec(3) = vz

  !! Initialise
  euler_rotate = 1;
  rvec(1) = 0;
  rvec(2) = 0;
  rvec(3) = 0;
 
  !! Create relevant rotation operator, Euler angles
  select case (shock_direction)
  case ("+x")
     !! Do nothing
     return
  case ("-x")
     !! alpha=PI, beta=gam=0
     rvec(1)= -ivec(1)
     rvec(2)= -ivec(2)
     rvec(3)=  ivec(3)
     euler_rotate=0;
  case ("+y")
     !! alpha=PI/2, beta=gam=0
     rvec(1)= -ivec(2)*vparity
     rvec(2)=  ivec(1)
     rvec(3)=  ivec(3)
     euler_rotate=0;
  case ("-y")
     !! alpha=(3/2)*PI, beta=gam=0
     rvec(1)= -ivec(2)
     rvec(2)= -ivec(1)
     rvec(3)=  ivec(3)
     euler_rotate=0;
  case ("+d")
     alpha=0.25d0*PI
     beta=0.d0
     gam=0.d0
  case ("-d")
     alpha=1.25d0*PI
     beta=0.d0
     gam=0.d0
  case ("+z")
     !! beta=-PI/2, alpha=gam=0
     rvec(1)= -ivec(3)
     rvec(2)=  ivec(2)
     rvec(3)=  ivec(1)
     euler_rotate=0;
  case ("-z")
     !! beta=3*PI/2, alpha=gam=0
     rvec(1)=  ivec(3)
     rvec(2)=  ivec(2)
     rvec(3)= -ivec(1)
     euler_rotate=0;
  case ("+r")
     if ( vx.ne.0d0 .or. vy.ne.0d0 .or. vz.ne.0d0 ) then
        call CCTK_WARN(0,"Please implement vector rotations for radial shock")
     end if
  case ("-r")
     if ( vx.ne.0d0 .or. vy.ne.0d0 .or. vz.ne.0d0 ) then
        call CCTK_WARN(0,"Please implement vector rotations for radial shock")
     end if
  case default
     call CCTK_WARN(0,"Unidentified shock case in Whisky_ID_Utils.")
  end select

  if ( euler_rotate.gt.0 ) then
     !! Create Euler rotation matrix
     Arot(1,1) = cos(alpha)*cos(beta)*cos(gam) - sin(alpha)*sin(gam)
     Arot(1,2) = -(cos(alpha)*cos(beta)*sin(gam)+sin(alpha)*cos(gam))
     Arot(1,3) = cos(alpha)*sin(beta)
     Arot(2,1) = sin(alpha)*cos(beta)*cos(gam) + cos(alpha)*sin(gam) 
     Arot(2,2) = -sin(alpha)*cos(beta)*sin(gam) + cos(alpha)*cos(gam)
     Arot(2,3) = sin(alpha)*sin(beta)
     Arot(3,1) = -sin(beta)*cos(gam)
     Arot(3,2) = sin(beta)*sin(gam)
     Arot(3,3) = cos(beta) 

     !! Rotate vector
     rvec(1)=0;
     rvec(2)=0;
     rvec(3)=0;
     do i=1,3
        do j=1,3
           rvec(i) = rvec(i) + Arot(i,j)*ivec(j) 
        end do 
     end do 
  end if

  !! Feed new vector back through original pointers
  vx=rvec(1) 
  vy=rvec(2) 
  vz=rvec(3) 

end subroutine RotateVector
