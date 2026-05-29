#include "cctk.h"
#include "cctk_Arguments.h"

      integer function Coordinates_Startup()

      implicit none

      integer ierr

c     Register coordinate systems

      print *
      print *,"CORRECTLY REGISTERING COORDINATE SYSTEM"
      print *,"---------------------------------------"

      print *
      print *,"CCTK_CoordRegisterSystem(ierr,3,<testcart3d>)"
      call CCTK_CoordRegisterSystem(ierr,3,"testcart3d")
      print *,"ierr=",ierr

      print *	
      print *,"CCTK_CoordRegisterSystem(ierr,3,<testspher3d>)"
      call CCTK_CoordRegisterSystem(ierr,3,"testspher3d")
      print *,"ierr=",ierr

      print *
      print *,"INCORRECTLY REGISTERING COORDINATE SYSTEM"
      print *,"-----------------------------------------"

      print *
      print *,"Registering duplicate system with same dimension..."
      print *,"CCTK_CoordRegisterSystem(ierr,3,<testcart3d>)"
      call CCTK_CoordRegisterSystem(ierr,3,"testcart3d")
      print *,"ierr=",ierr
	
      print *
      print *,"Registering duplicate system with different dimension..."
      print *,"CCTK_CoordRegisterSystem(ierr,2,<testcart3d>)"
      call CCTK_CoordRegisterSystem(ierr,2,"testcart3d")
      print *,"ierr=",ierr
	
      print *
      print *,"Registering with wrong dimension ..."
      print *,"CCTK_CoordRegisterSystem(ierr,0,<newone>)"
      call CCTK_CoordRegisterSystem(ierr,0,"newone")
      print *,"ierr=",ierr
	
      print *
      print *,"Registering with empty name ..."
      print *,"CCTK_CoordRegisterSystem(ierr,0,<>)"
      call CCTK_CoordRegisterSystem(ierr,0," ")
      print *,"ierr=",ierr
	
      Coordinates_Startup = 0
	
      return 
      end



      subroutine Coordinates_Setup(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer i,j,k,ierr
      CCTK_REAL lower,upper

      do k=cctk_lbnd(3),cctk_ubnd(3)
	 do j=cctk_lbnd(2),cctk_ubnd(2)
	    do i=cctk_lbnd(1),cctk_ubnd(1)
	       x(i,j,k) = i
	       y(i,j,k) = j
	       z(i,j,k) = k
            end do
         end do
      end do

      do k=cctk_lbnd(3),cctk_ubnd(3)
	 do j=cctk_lbnd(2),cctk_ubnd(2)
	    do i=cctk_lbnd(1),cctk_ubnd(1)
	       r(i,j,k) = sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
            end do
         end do
      end do

      print *
      print *,"REGISTERING COORDINATES"
      print *,"-----------------------"
      print *,"CCTK_CoordRegisterData(ierr,1,<coordinates::x>,<x>,<testcart3d>)"
      call CCTK_CoordRegisterData(ierr,1,"coordinates::x","x","testcart3d")
      print *,"ierr = ",ierr
      print *,"CCTK_CoordRegisterData(ierr,2,<coordinates::y>,<y>,<testcart3d>)"
      call CCTK_CoordRegisterData(ierr,2,"coordinates::y","y","testcart3d")
      print *,"ierr = ",ierr

      print *
      print *,"WRONGLY REGISTERING A COORDINATE"
      print *,"--------------------------------"
	
      print *,"Registering for a non-existent system ..."
      print *,"CCTK_CoordRegisterData(ierr,1,<coordinates::x>,<new>,<silly>)"
      call CCTK_CoordRegisterData(ierr,1,"coordinates::x","new","silly")
      print *,"ierr = ",ierr

      print *,"Registering outside dimension ..."
      print *,"CCTK_CoordRegisterData(ierr,4,<coordinates::x>,<new>,<testcart3d>)"
      call CCTK_CoordRegisterData(ierr,4,"coordinates::x","new","testcart3d")
      print *,"ierr = ",ierr

      print *,"Registering for the same direction ..."
      print *,"CCTK_CoordRegisterData(ierr,1,<coordinates::x>,<new>,<testcart3d>)"
      call CCTK_CoordRegisterData(ierr,1,"coordinates::x","new","testcart3d")
      print *,"ierr = ",ierr

      print *,"Registering with the same name ..."
      print *,"CCTK_CoordRegisterData(ierr,3,<coordinates::x>,<x>,<testcart3d>)"
      call CCTK_CoordRegisterData(ierr,3,"coordinates::x","x","testcart3d")
      print *,"ierr = ",ierr

      print *,"Registering something not a grid variable ..."
      print *,"CCTK_CoordRegisterData(ierr,3,<coordinates::silly>,<z>,<testcart3d>)"
      call CCTK_CoordRegisterData(ierr,3,"coordinates::silly","z","testcart3d")
      print *,"ierr = ",ierr


      print *
      print *,"CORRECTLY REGISTERING COORDINATE RANGE"
      print *,"--------------------------------------"

c     Register range for x
      lower = 1.0
      upper = cctk_gsh(1)
      print *
      print *,"CCTK_CoordRegisterRange(ierr,cctkGH,lower,upper,1,<>,<testcart3d>)"
      call CCTK_CoordRegisterRange(ierr,cctkGH,lower,upper,1,"","testcart3d")
      print *,"ierr=",ierr

c     Register range for y
      lower = 1.0
      upper = cctk_gsh(2)
      print *
      print *,"CCTK_CoordRegisterRange(ierr,cctkGH,lower,upper,-1,<y>,<testcart3d>)"
      call CCTK_CoordRegisterRange(ierr,cctkGH,lower,upper,-1,"y","testcart3d")
      print *,"ierr=",ierr

      print *
      print *,"WRONGLY REGISTERING COORDINATE RANGE"
      print *,"------------------------------------"

c     Registering another range for the same direction
      print *,"Registering with the same direction ..."
      lower = -1.0
      upper = 1.0
      print *,"CCTK_CoordRegisterRange(ierr,cctkGH,lower,upper,1,<>,<testcart3d>)"
      call CCTK_CoordRegisterRange(ierr,cctkGH,lower,upper,1,"","testcart3d")
      print *,"ierr=",ierr

      end



      subroutine Coordinates_Query(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer dim,index,ierr,dir
      CCTK_REAL lower,upper

      print *
      print *,"WRONGLY QUERYING DIMENSION"
      print *,"----------------------------"

      print *
      print *,"Wrong system name"
      print *,"CCTK_CoordSystemDim(dim,<silly>)"
      call CCTK_CoordSystemDim(dim,"silly")
      print *,"dim=",dim

      print *
      print *,"CORRECTLY QUERYING DIMENSION"
      print *,"----------------------------"

      print *
      print *,"CCTK_CoordSystemDim(dim,<testcart3d>)"
      call CCTK_CoordSystemDim(dim,"testcart3d")
      print *,"dim=",dim

      print *
      print *,"CCTK_CoordSystemDim(dim,<testspher3d>)"
      call CCTK_CoordSystemDim(dim,"testspher3d")
      print *,"dim=",dim

      print *
      print *,"CORRECTLY QUERYING INDEX"
      print *,"------------------------"

      print *
      print *,"CCTK_CoordIndex(index,1,<>,<testcart3d>)"
      call CCTK_CoordIndex(index,1,"","testcart3d")
      print *,"index=",index

      print *
      print *,"CCTK_CoordIndex(index,2,<>,<testcart3d>)"
      call CCTK_CoordIndex(index,2,"","testcart3d")
      print *,"index=",index

      print *
      print *,"CCTK_CoordIndex(index,-1,<x>,<testcart3d>)"
      call CCTK_CoordIndex(index,-1,"x","testcart3d")
      print *,"index=",index

      print *
      print *,"CCTK_CoordIndex(index,-1,<y>,<testcart3d>)"
      call CCTK_CoordIndex(index,-1,"y","testcart3d")
      print *,"index=",index

      print *
      print *,"WRONGLY QUERYING INDEX"
      print *,"----------------------"
     
      print *
      print *,"Wrong system name ..."
      print *,"CCTK_CoordIndex(index,1,<>,<silly>)"
      call CCTK_CoordIndex(index,1,"","silly")
      print *,"index=",index
      
      print *
      print *,"Direction out of range ..."
      print *,"CCTK_CoordIndex(index,4,<>,<testcart3d>)"
      call CCTK_CoordIndex(index,4,"","testcart3d")
      print *,"index=",index

      print *
      print *,"Wrong coordinate name ..."
      print *,"CCTK_CoordIndex(index,-1,<silly>,<testcart3d>)"
      call CCTK_CoordIndex(index,-1,"silly","testcart3d")
      print *,"index=",index

      print *
      print *,"Unregistered coordinate ..."
      print *,"CCTK_CoordIndex(index,3,<>,<testcart3d>)"
      call CCTK_CoordIndex(index,3,"","testcart3d")
      print *,"index=",index


      print *
      print *,"WRONGLY QUERYING RANGE"
      print *,"----------------------"

      print *,"Wrong system name ..."
      print *,"CCTK_CoordRange(ierr,cctkGH,lower,upper,1,<>,<silly>)"
      call CCTK_CoordRange(ierr,cctkGH,lower,upper,1,"","silly")
      print *,"ierr=",ierr
      print *,"lower=",lower
      print *,"upper=",upper

      print *
      print *,"Wrong direction ..."
      print *,"CCTK_CoordRange(ierr,cctkGH,lower,upper,34,<>,<cartgrid3d>)"
      call CCTK_CoordRange(ierr,cctkGH,lower,upper,1,"","cartgrid3d")
      print *,"ierr=",ierr
      print *,"lower=",lower
      print *,"upper=",upper

      print *
      print *,"Range not registered ..."
      call CCTK_CoordRange(ierr,cctkGH,lower,upper,3,"","cartgrid3d")
      print *,"CCTK_CoordRange(ierr,cctkGH,lower,upper,3,<>,<cartgrid3d>)"
      print *,"ierr=",ierr
      print *,"lower=",lower
      print *,"upper=",upper


      print *
      print *,"CORRECTLY QUERYING RANGE"
      print *,"------------------------"

      print *
      print *,"CCTK_CoordRange(ierr,cctkGH,lower,upper,1,<>,<testcart3d>)"
      call CCTK_CoordRange(ierr,cctkGH,lower,upper,1,"","testcart3d")
      print *,"ierr=",ierr
      print *,"lower=",lower
      print *,"upper=",upper

      print *
      print *,"CCTK_CoordRange(ierr,cctkGH,lower,upper,2,<>,<testcart3d>)"
      call CCTK_CoordRange(ierr,cctkGH,lower,upper,2,"","testcart3d")
      print *,"ierr=",ierr
      print *,"lower=",lower
      print *,"upper=",upper


      print *
      print *,"INCORRECTLY QUERYING COORDINATE DIRECTION"
      print *,"-----------------------------------------"

      print *
      print *,"Wrong system name ..."
      print *,"CCTK_CoordDir(dir,<x>,<silly>)"
      call CCTK_CoordDir(dir,"x","silly")
      print *,"dir=",dir

      print *
      print *,"Wrong coordinate name ..."
      print *,"CCTK_CoordDir(dir,<silly>,<testcart3d>)"
      call CCTK_CoordDir(dir,"silly","testcart3d")
      print *,"dir=",dir

      print *
      print *,"CORRECTLY QUERYING COORDINATE DIRECTION"
      print *,"---------------------------------------"

      print *
      print *,"CCTK_CoordDir(dir,<x>,<testcart3d>)"
      call CCTK_CoordDir(dir,"x","testcart3d")
      print *,"dir=",dir

      end
