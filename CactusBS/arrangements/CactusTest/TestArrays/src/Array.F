 /*@@
   @file      Array.F
   @date      
   @author    Gabrielle Allen
   @desc 
              Set up and use arrays
   @enddesc 
 @@*/
  

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine SetUpFortranArrays(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k,m,ierr
      integer myproc
      integer lsh_1d(1),lsh_2d(2),lsh_3d(3),lsh_4d(4)
      integer gsh_1d(1),gsh_2d(2),gsh_3d(3),gsh_4d(4)
      integer lbnd_1d(1),lbnd_2d(2),lbnd_3d(3),lbnd_4d(4)
      integer ubnd_1d(1),ubnd_2d(2),ubnd_3d(3),ubnd_4d(4)
      integer nghostzones_1d(1),nghostzones_2d(2),nghostzones_3d(3),nghostzones_4d(4)
      integer bbox_1d(2),bbox_2d(4),bbox_3d(6),bbox_4d(8)
      integer varindex1,varindex2,varindex3,varindex4
      integer groupindex1,groupindex2,groupindex3,groupindex4

      print *
      print *,"In SetUpFortranArrays"
      print *,"---------------------"
      print *

      call CCTK_VarIndex(varindex1,"testarrays::gf1d")
      call CCTK_VarIndex(varindex2,"testarrays::gf2d")
      call CCTK_VarIndex(varindex3,"testarrays::gf3d")
      call CCTK_VarIndex(varindex4,"testarrays::gf4d")
      call CCTK_GroupIndex(groupindex1,"testarrays::foo1d")
      call CCTK_GroupIndex(groupindex2,"testarrays::foo2d")
      call CCTK_GroupIndex(groupindex3,"testarrays::foo3d")
      call CCTK_GroupIndex(groupindex4,"testarrays::foo4d")

      print *,"Checking array queries"
      print *,"----------------------"

      print *
      print *,"Local shape (GN/VN/GI/VI) ..."
      call CCTK_GrouplshGN(ierr,cctkGH,1,lsh_1d,"testarrays::foo1d")
      print *,lsh_1d(1)
      call CCTK_GrouplshVN(ierr,cctkGH,1,lsh_1d,"testarrays::gf1d")
      print *,lsh_1d(1)
      call CCTK_GrouplshGI(ierr,cctkGH,1,lsh_1d,groupindex1)
      print *,lsh_1d(1)
      call CCTK_GrouplshVI(ierr,cctkGH,1,lsh_1d,varindex1)
      print *,lsh_1d(1)

      call CCTK_GrouplshGN(ierr,cctkGH,2,lsh_2d,"testarrays::foo2d")
      print *,lsh_2d(1),lsh_2d(2)
      call CCTK_GrouplshVN(ierr,cctkGH,2,lsh_2d,"testarrays::gf2d")
      print *,lsh_2d(1),lsh_2d(2)
      call CCTK_GrouplshGI(ierr,cctkGH,2,lsh_2d,groupindex2)
      print *,lsh_2d(1),lsh_2d(2)
      call CCTK_GrouplshVI(ierr,cctkGH,2,lsh_2d,varindex2)
      print *,lsh_2d(1),lsh_2d(2)

      call CCTK_GrouplshGN(ierr,cctkGH,3,lsh_3d,"testarrays::foo3d")
      print *,lsh_3d(1),lsh_3d(2),lsh_3d(3)
      call CCTK_GrouplshVN(ierr,cctkGH,3,lsh_3d,"testarrays::gf3d")
      print *,lsh_3d(1),lsh_3d(2),lsh_3d(3)
      call CCTK_GrouplshGI(ierr,cctkGH,3,lsh_3d,groupindex3)
      print *,lsh_3d(1),lsh_3d(2),lsh_3d(3)
      call CCTK_GrouplshVI(ierr,cctkGH,3,lsh_3d,varindex3)
      print *,lsh_3d(1),lsh_3d(2),lsh_3d(3)

      call CCTK_GrouplshGN(ierr,cctkGH,4,lsh_4d,"testarrays::foo4d")
      print *,lsh_4d(1),lsh_4d(2),lsh_4d(3),lsh_4d(4)
      call CCTK_GrouplshVN(ierr,cctkGH,4,lsh_4d,"testarrays::gf4d")
      print *,lsh_4d(1),lsh_4d(2),lsh_4d(3),lsh_4d(4)
      call CCTK_GrouplshGI(ierr,cctkGH,4,lsh_4d,groupindex4)
      print *,lsh_4d(1),lsh_4d(2),lsh_4d(3),lsh_4d(4)
      call CCTK_GrouplshVI(ierr,cctkGH,4,lsh_4d,varindex4)
      print *,lsh_4d(1),lsh_4d(2),lsh_4d(3),lsh_4d(4)

      print *
      print *,"Global shape (GN/VN/GI/VI) ..."
      call CCTK_GroupgshGN(ierr,cctkGH,1,gsh_1d,"testarrays::foo1d")
      print *,gsh_1d(1)
      call CCTK_GroupgshVN(ierr,cctkGH,1,gsh_1d,"testarrays::gf1d")
      print *,gsh_1d(1)
      call CCTK_GroupgshGI(ierr,cctkGH,1,gsh_1d,groupindex1)
      print *,gsh_1d(1)
      call CCTK_GroupgshVI(ierr,cctkGH,1,gsh_1d,varindex1)
      print *,gsh_1d(1)

      call CCTK_GroupgshGN(ierr,cctkGH,2,gsh_2d,"testarrays::foo2d")
      print *,gsh_2d(1),gsh_2d(2)
      call CCTK_GroupgshVN(ierr,cctkGH,2,gsh_2d,"testarrays::gf2d")
      print *,gsh_2d(1),gsh_2d(2)
      call CCTK_GroupgshGI(ierr,cctkGH,2,gsh_2d,groupindex2)
      print *,gsh_2d(1),gsh_2d(2)
      call CCTK_GroupgshVI(ierr,cctkGH,2,gsh_2d,varindex2)
      print *,gsh_2d(1),gsh_2d(2)

      call CCTK_GroupgshGN(ierr,cctkGH,3,gsh_3d,"testarrays::foo3d")
      print *,gsh_3d(1),gsh_3d(2),gsh_3d(3)
      call CCTK_GroupgshVN(ierr,cctkGH,3,gsh_3d,"testarrays::gf3d")
      print *,gsh_3d(1),gsh_3d(2),gsh_3d(3)
      call CCTK_GroupgshGI(ierr,cctkGH,3,gsh_3d,groupindex3)
      print *,gsh_3d(1),gsh_3d(2),gsh_3d(3)
      call CCTK_GroupgshVI(ierr,cctkGH,3,gsh_3d,varindex3)
      print *,gsh_3d(1),gsh_3d(2),gsh_3d(3)

      call CCTK_GroupgshGN(ierr,cctkGH,4,gsh_4d,"testarrays::foo4d")
      print *,gsh_4d(1),gsh_4d(2),gsh_4d(3),gsh_4d(4)
      call CCTK_GroupgshVN(ierr,cctkGH,4,gsh_4d,"testarrays::gf4d")
      print *,gsh_4d(1),gsh_4d(2),gsh_4d(3),gsh_4d(4)
      call CCTK_GroupgshGI(ierr,cctkGH,4,gsh_4d,groupindex4)
      print *,gsh_4d(1),gsh_4d(2),gsh_4d(3),gsh_4d(4)
      call CCTK_GroupgshVI(ierr,cctkGH,4,gsh_4d,varindex4)
      print *,gsh_4d(1),gsh_4d(2),gsh_4d(3),gsh_4d(4)

      print *
      print *,"Lower bounds (GN/VN/GI/VI) ..."
      call CCTK_GrouplbndGN(ierr,cctkGH,1,lbnd_1d,"testarrays::foo1d")
      print *,lbnd_1d(1)
      call CCTK_GrouplbndVN(ierr,cctkGH,1,lbnd_1d,"testarrays::gf1d")
      print *,lbnd_1d(1)
      call CCTK_GrouplbndGI(ierr,cctkGH,1,lbnd_1d,groupindex1)
      print *,lbnd_1d(1)
      call CCTK_GrouplbndVI(ierr,cctkGH,1,lbnd_1d,varindex1)
      print *,lbnd_1d(1)

      call CCTK_GrouplbndGN(ierr,cctkGH,2,lbnd_2d,"testarrays::foo2d")
      print *,lbnd_2d(1),lbnd_2d(2)
      call CCTK_GrouplbndVN(ierr,cctkGH,2,lbnd_2d,"testarrays::gf2d")
      print *,lbnd_2d(1),lbnd_2d(2)
      call CCTK_GrouplbndGI(ierr,cctkGH,2,lbnd_2d,groupindex2)
      print *,lbnd_2d(1),lbnd_2d(2)
      call CCTK_GrouplbndVI(ierr,cctkGH,2,lbnd_2d,varindex2)
      print *,lbnd_2d(1),lbnd_2d(2)

      call CCTK_GrouplbndGN(ierr,cctkGH,3,lbnd_3d,"testarrays::foo3d")
      print *,lbnd_3d(1),lbnd_3d(2),lbnd_3d(3)
      call CCTK_GrouplbndVN(ierr,cctkGH,3,lbnd_3d,"testarrays::gf3d")
      print *,lbnd_3d(1),lbnd_3d(2),lbnd_3d(3)
      call CCTK_GrouplbndGI(ierr,cctkGH,3,lbnd_3d,groupindex3)
      print *,lbnd_3d(1),lbnd_3d(2),lbnd_3d(3)
      call CCTK_GrouplbndVI(ierr,cctkGH,3,lbnd_3d,varindex3)
      print *,lbnd_3d(1),lbnd_3d(2),lbnd_3d(3)

      call CCTK_GrouplbndGN(ierr,cctkGH,4,lbnd_4d,"testarrays::foo4d")
      print *,lbnd_4d(1),lbnd_4d(2),lbnd_4d(3),lbnd_4d(4)
      call CCTK_GrouplbndVN(ierr,cctkGH,4,lbnd_4d,"testarrays::gf4d")
      print *,lbnd_4d(1),lbnd_4d(2),lbnd_4d(3),lbnd_4d(4)
      call CCTK_GrouplbndGI(ierr,cctkGH,4,lbnd_4d,groupindex4)
      print *,lbnd_4d(1),lbnd_4d(2),lbnd_4d(3),lbnd_4d(4)
      call CCTK_GrouplbndVI(ierr,cctkGH,4,lbnd_4d,varindex4)
      print *,lbnd_4d(1),lbnd_4d(2),lbnd_4d(3),lbnd_4d(4)

      print *
      print *,"Upper bounds (GN/VN/GI/VI) ..."
      call CCTK_GroupubndGN(ierr,cctkGH,1,ubnd_1d,"testarrays::foo1d")
      print *,ubnd_1d(1)
      call CCTK_GroupubndVN(ierr,cctkGH,1,ubnd_1d,"testarrays::gf1d")
      print *,ubnd_1d(1)
      call CCTK_GroupubndGI(ierr,cctkGH,1,ubnd_1d,groupindex1)
      print *,ubnd_1d(1)
      call CCTK_GroupubndVI(ierr,cctkGH,1,ubnd_1d,varindex1)
      print *,ubnd_1d(1)

      call CCTK_GroupubndGN(ierr,cctkGH,2,ubnd_2d,"testarrays::foo2d")
      print *,ubnd_2d(1),ubnd_2d(2)
      call CCTK_GroupubndVN(ierr,cctkGH,2,ubnd_2d,"testarrays::gf2d")
      print *,ubnd_2d(1),ubnd_2d(2)
      call CCTK_GroupubndGI(ierr,cctkGH,2,ubnd_2d,groupindex2)
      print *,ubnd_2d(1),ubnd_2d(2)
      call CCTK_GroupubndVI(ierr,cctkGH,2,ubnd_2d,varindex2)
      print *,ubnd_2d(1),ubnd_2d(2)

      call CCTK_GroupubndGN(ierr,cctkGH,3,ubnd_3d,"testarrays::foo3d")
      print *,ubnd_3d(1),ubnd_3d(2),ubnd_3d(3)
      call CCTK_GroupubndVN(ierr,cctkGH,3,ubnd_3d,"testarrays::gf3d")
      print *,ubnd_3d(1),ubnd_3d(2),ubnd_3d(3)
      call CCTK_GroupubndGI(ierr,cctkGH,3,ubnd_3d,groupindex3)
      print *,ubnd_3d(1),ubnd_3d(2),ubnd_3d(3)
      call CCTK_GroupubndVI(ierr,cctkGH,3,ubnd_3d,varindex3)
      print *,ubnd_3d(1),ubnd_3d(2),ubnd_3d(3)

      call CCTK_GroupubndGN(ierr,cctkGH,4,ubnd_4d,"testarrays::foo4d")
      print *,ubnd_4d(1),ubnd_4d(2),ubnd_4d(3),ubnd_4d(4)
      call CCTK_GroupubndVN(ierr,cctkGH,4,ubnd_4d,"testarrays::gf4d")
      print *,ubnd_4d(1),ubnd_4d(2),ubnd_4d(3),ubnd_4d(4)
      call CCTK_GroupubndGI(ierr,cctkGH,4,ubnd_4d,groupindex4)
      print *,ubnd_4d(1),ubnd_4d(2),ubnd_4d(3),ubnd_4d(4)
      call CCTK_GroupubndVI(ierr,cctkGH,4,ubnd_4d,varindex4)
      print *,ubnd_4d(1),ubnd_4d(2),ubnd_4d(3),ubnd_4d(4)

      print *
      print *,"Ghostzones (GN/VN/GI/VI) ..."
      call CCTK_GroupnghostzonesGN(ierr,cctkGH,1,nghostzones_1d,"testarrays::foo1d")
      print *,nghostzones_1d(1)
      call CCTK_GroupnghostzonesVN(ierr,cctkGH,1,nghostzones_1d,"testarrays::gf1d")
      print *,nghostzones_1d(1)
      call CCTK_GroupnghostzonesGI(ierr,cctkGH,1,nghostzones_1d,groupindex1)
      print *,nghostzones_1d(1)
      call CCTK_GroupnghostzonesVI(ierr,cctkGH,1,nghostzones_1d,varindex1)
      print *,nghostzones_1d(1)

      call CCTK_GroupnghostzonesGN(ierr,cctkGH,2,nghostzones_2d,"testarrays::foo2d")
      print *,nghostzones_2d(1),nghostzones_2d(2)
      call CCTK_GroupnghostzonesVN(ierr,cctkGH,2,nghostzones_2d,"testarrays::gf2d")
      print *,nghostzones_2d(1),nghostzones_2d(2)
      call CCTK_GroupnghostzonesGI(ierr,cctkGH,2,nghostzones_2d,groupindex2)
      print *,nghostzones_2d(1),nghostzones_2d(2)
      call CCTK_GroupnghostzonesVI(ierr,cctkGH,2,nghostzones_2d,varindex2)
      print *,nghostzones_2d(1),nghostzones_2d(2)

      call CCTK_GroupnghostzonesGN(ierr,cctkGH,3,nghostzones_3d,"testarrays::foo3d")
      print *,nghostzones_3d(1),nghostzones_3d(2),nghostzones_3d(3)
      call CCTK_GroupnghostzonesVN(ierr,cctkGH,3,nghostzones_3d,"testarrays::gf3d")
      print *,nghostzones_3d(1),nghostzones_3d(2),nghostzones_3d(3)
      call CCTK_GroupnghostzonesGI(ierr,cctkGH,3,nghostzones_3d,groupindex3)
      print *,nghostzones_3d(1),nghostzones_3d(2),nghostzones_3d(3)
      call CCTK_GroupnghostzonesVI(ierr,cctkGH,3,nghostzones_3d,varindex3)
      print *,nghostzones_3d(1),nghostzones_3d(2),nghostzones_3d(3)

      call CCTK_GroupnghostzonesGN(ierr,cctkGH,4,nghostzones_4d,"testarrays::foo4d")
      print *,nghostzones_4d(1),nghostzones_4d(2),nghostzones_4d(3),nghostzones_4d(4)
      call CCTK_GroupnghostzonesVN(ierr,cctkGH,4,nghostzones_4d,"testarrays::gf4d")
      print *,nghostzones_4d(1),nghostzones_4d(2),nghostzones_4d(3),nghostzones_4d(4)
      call CCTK_GroupnghostzonesGI(ierr,cctkGH,4,nghostzones_4d,groupindex4)
      print *,nghostzones_4d(1),nghostzones_4d(2),nghostzones_4d(3),nghostzones_4d(4)
      call CCTK_GroupnghostzonesVI(ierr,cctkGH,4,nghostzones_4d,varindex4)
      print *,nghostzones_4d(1),nghostzones_4d(2),nghostzones_4d(3),nghostzones_4d(4)

      print *
      print *,"Bounding box (GN/VN/GI/VI) ..."
      call CCTK_GroupbboxGN(ierr,cctkGH,2,bbox_1d,"testarrays::foo1d")
      print *,bbox_1d(1),bbox_1d(2)
      call CCTK_GroupbboxVN(ierr,cctkGH,2,bbox_1d,"testarrays::gf1d")
      print *,bbox_1d(1),bbox_1d(2)
      call CCTK_GroupbboxGI(ierr,cctkGH,2,bbox_1d,groupindex1)
      print *,bbox_1d(1),bbox_1d(2)
      call CCTK_GroupbboxVI(ierr,cctkGH,2,bbox_1d,varindex1)
      print *,bbox_1d(1),bbox_1d(2)

      call CCTK_GroupbboxGN(ierr,cctkGH,4,bbox_2d,"testarrays::foo2d")
      print *,bbox_2d(1),bbox_2d(2),bbox_2d(3),bbox_2d(4)
      call CCTK_GroupbboxVN(ierr,cctkGH,4,bbox_2d,"testarrays::gf2d")
      print *,bbox_2d(1),bbox_2d(2),bbox_2d(3),bbox_2d(4)
      call CCTK_GroupbboxGI(ierr,cctkGH,4,bbox_2d,groupindex2)
      print *,bbox_2d(1),bbox_2d(2),bbox_2d(3),bbox_2d(4)
      call CCTK_GroupbboxVI(ierr,cctkGH,4,bbox_2d,varindex2)
      print *,bbox_2d(1),bbox_2d(2),bbox_2d(3),bbox_2d(4)

      call CCTK_GroupbboxGN(ierr,cctkGH,6,bbox_3d,"testarrays::foo3d")
      print *,bbox_3d(1),bbox_3d(2),bbox_3d(3),
     &     bbox_3d(4),bbox_3d(5),bbox_3d(6)
      call CCTK_GroupbboxVN(ierr,cctkGH,6,bbox_3d,"testarrays::gf3d")
      print *,bbox_3d(1),bbox_3d(2),bbox_3d(3),
     &     bbox_3d(4),bbox_3d(5),bbox_3d(6)
      call CCTK_GroupbboxGI(ierr,cctkGH,6,bbox_3d,groupindex3)
      print *,bbox_3d(1),bbox_3d(2),bbox_3d(3),
     &     bbox_3d(4),bbox_3d(5),bbox_3d(6)
      call CCTK_GroupbboxVI(ierr,cctkGH,6,bbox_3d,varindex3)
      print *,bbox_3d(1),bbox_3d(2),bbox_3d(3),
     &     bbox_3d(4),bbox_3d(5),bbox_3d(6)

      call CCTK_GroupbboxGN(ierr,cctkGH,8,bbox_4d,"testarrays::foo4d")
      print *,bbox_4d(1),bbox_4d(2),bbox_4d(3),
     &     bbox_4d(4),bbox_4d(5),bbox_4d(6),bbox_4d(7),bbox_4d(8)
      call CCTK_GroupbboxVN(ierr,cctkGH,8,bbox_4d,"testarrays::gf4d")
      print *,bbox_4d(1),bbox_4d(2),bbox_4d(3),
     &     bbox_4d(4),bbox_4d(5),bbox_4d(6),bbox_4d(7),bbox_4d(8)
      call CCTK_GroupbboxGI(ierr,cctkGH,8,bbox_4d,groupindex4)
      print *,bbox_4d(1),bbox_4d(2),bbox_4d(3),
     &     bbox_4d(4),bbox_4d(5),bbox_4d(6),bbox_4d(7),bbox_4d(8)
      call CCTK_GroupbboxVI(ierr,cctkGH,8,bbox_4d,varindex4)
      print *,bbox_4d(1),bbox_4d(2),bbox_4d(3),
     &     bbox_4d(4),bbox_4d(5),bbox_4d(6),bbox_4d(7),bbox_4d(8)

      print *
      print *,"1D Array size is ",lsh_1d(1)
      do i=1,lsh_1d(1)
         gf1d(i) = dble(lbnd_1d(1)+i)
      end do

      print *,"2D Array size is ",lsh_2d(1),lsh_2d(2)
      do i=1,lsh_2d(1)
         do j=1,lsh_2d(2)
               gf2d(i,j) = (lbnd_2d(1)+i)*(lbnd_2d(2)+j)
         end do
      end do

      print *,"3D Array size is ",lsh_3d(1),lsh_3d(2),lsh_3d(3)
      do i=1,lsh_3d(1)
         do j=1,lsh_3d(2)
            do k=1,lsh_3d(3)
               gf3d(i,j,k) = (lbnd_3d(1)+i)*(lbnd_3d(2)+j)*
     &              (lbnd_3d(3)+k)
            end do
         end do
      end do

      print *,"4D Array size is ",lsh_4d(1),lsh_4d(2),lsh_4d(3),lsh_4d(4)
      do i=1,lsh_4d(1)
         do j=1,lsh_4d(2)
            do k=1,lsh_4d(3)
               do m=1,lsh_4d(4)
                 gf4d(i,j,k,m) = (lbnd_4d(1)+i)*(lbnd_4d(2)+j)*
     &                (lbnd_4d(3)+k)*(lbnd_4d(4)+m)
               end do
            end do
         end do
      end do
      
      end subroutine SetUpFortranArrays


      subroutine UseFortranArrays(CCTK_ARGUMENTS)
      
      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k,m,ierr
      integer lsh_1d(1),lsh_2d(2),lsh_3d(3),lsh_4d(4)
      integer nghostzones_1d(1),nghostzones_2d(2),nghostzones_3d(3),nghostzones_4d(4)

      print *
      print *,"In UseFortranArrays"
      print *,"-------------------"

      call CCTK_GrouplshGN(ierr,cctkGH,1,lsh_1d,"testarrays::foo1d")
      call CCTK_GrouplshGN(ierr,cctkGH,2,lsh_2d,"testarrays::foo2d")
      call CCTK_GrouplshGN(ierr,cctkGH,3,lsh_3d,"testarrays::foo3d")
      call CCTK_GrouplshGN(ierr,cctkGH,4,lsh_4d,"testarrays::foo4d")

      call CCTK_GroupnghostzonesGN(ierr,cctkGH,1,nghostzones_1d,"testarrays::foo1d")
      call CCTK_GroupnghostzonesGN(ierr,cctkGH,2,nghostzones_2d,"testarrays::foo2d")
      call CCTK_GroupnghostzonesGN(ierr,cctkGH,3,nghostzones_3d,"testarrays::foo3d")
      call CCTK_GroupnghostzonesGN(ierr,cctkGH,4,nghostzones_4d,"testarrays::foo4d")

c     1D Array
c     ========
      print *,"1D Array size is ",lsh_1d(1)
      do i=1+nghostzones_1d(1),lsh_1d(1)-nghostzones_1d(1)
         gf1d(i)=gf1d(i)+1
      end do
      print *,"Syncing 1D array ..."
      call CCTK_SyncGroup(ierr,cctkGH,"testarrays::foo1d")

c     2D Array
c     ========
      print *,"2D Array size is ",lsh_2d(1),lsh_2d(2)
      do i=1+nghostzones_2d(1),lsh_2d(1)-nghostzones_2d(1)
         do j=1+nghostzones_2d(2),lsh_2d(2)-nghostzones_2d(2)
            gf2d(i,j)=gf2d(i,j)+1
         end do
      end do
      print *,"Syncing 2D array ..."
      call CCTK_SyncGroup(ierr,cctkGH,"testarrays::foo2d")
      
c     3D ARRAY
c     ========
      print *,"3D Array size is ",lsh_3d(1),lsh_3d(2),lsh_3d(3)
      do i=1+nghostzones_3d(1),lsh_3d(1)-nghostzones_3d(1)
         do j=1+nghostzones_3d(2),lsh_3d(2)-nghostzones_3d(2)
            do k=1+nghostzones_3d(3),lsh_3d(3)-nghostzones_3d(3)
               gf3d(i,j,k)=gf3d(i,j,k)+1
            end do
         end do
      end do
      print *,"Syncing 3D array ..."
      call CCTK_SyncGroup(ierr,cctkGH,"testarrays::foo3d")
      print *
      
c     4D ARRAY
c     ========
      print *,"4D Array size is ",lsh_4d(1),lsh_4d(2),lsh_4d(3),lsh_4d(4)
      do i=1+nghostzones_4d(1),lsh_4d(1)-nghostzones_4d(1)
         do j=1+nghostzones_4d(2),lsh_4d(2)-nghostzones_4d(2)
            do k=1+nghostzones_4d(3),lsh_4d(3)-nghostzones_4d(3)
               do m=1+nghostzones_4d(4),lsh_4d(4)-nghostzones_4d(4)
                  gf4d(i,j,k,m)=gf4d(i,j,k,m)+1
               end do
            end do
         end do
      end do
      print *,"Syncing 4D array ..."
      call CCTK_SyncGroup(ierr,cctkGH,"testarrays::foo4d")
      print *

      end subroutine UseFortranArrays

