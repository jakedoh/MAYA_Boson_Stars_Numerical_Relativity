! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine TestFreeF90 (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer a

!!$  ! Ensure that trailing commas are removed
!!$  call dummy (,)
!!$  call dummy (, )
!!$  call dummy (,,)
!!$  call dummy (,, )
!!$  call dummy (, ,)
!!$  call dummy (, , )
!!$  call dummy (,,,)
!!$  call dummy (,,, )
!!$  call dummy (,, ,)
!!$  call dummy (,, , )
!!$  call dummy (, ,,)
!!$  call dummy (, ,, )
!!$  call dummy (, , ,)
!!$  call dummy (, , , )
  
  ! Lines with varying length
  a =                                                                    75
  a =                                                                     76
  a =                                                                      77
  a =                                                                       78
  a =                                                                        79
  a =                                                                         80
  a =                                                                          81
  a =                                                                           82
  a =                                                                            83
  a =                                                                             84
  a =                                                                              85
  
  ! Lines with varying lengths and with a trailing space
  a =                                                                   75 
  a =                                                                    76 
  a =                                                                     77 
  a =                                                                      78 
  a =                                                                       79 
  a =                                                                        80 
  a =                                                                         81 
  a =                                                                          82 
  a =                                                                           83 
  a =                                                                            84 
  a =                                                                             85 
  
  ! Lines with varying lengths and with two trailing spaces
  a =                                                                  75  
  a =                                                                   76  
  a =                                                                    77  
  a =                                                                     78  
  a =                                                                      79  
  a =                                                                       80  
  a =                                                                        81  
  a =                                                                         82  
  a =                                                                          83  
  a =                                                                           84  
  a =                                                                            85  
  
  ! Lines with varying lengths and ending in a comment
  a =                                                                   75!
  a =                                                                    76!
  a =                                                                     77!
  a =                                                                      78!
  a =                                                                       79!
  a =                                                                        80!
  a =                                                                         81!
  a =                                                                          82!
  a =                                                                           83!
  a =                                                                            84!
  a =                                                                             85!
  
  ! Lines with varying lengths and ending in a comment plus a space
  a =                                                                  75! 
  a =                                                                   76! 
  a =                                                                    77! 
  a =                                                                     78! 
  a =                                                                      79! 
  a =                                                                       80! 
  a =                                                                        81! 
  a =                                                                         82! 
  a =                                                                          83! 
  a =                                                                           84! 
  a =                                                                            85! 
  
  ! Lines with varying lengths and ending in a comment plus two spaces
  a =                                                                 75!  
  a =                                                                  76!  
  a =                                                                   77!  
  a =                                                                    78!  
  a =                                                                     79!  
  a =                                                                      80!  
  a =                                                                       81!  
  a =                                                                        82!  
  a =                                                                         83!  
  a =                                                                          84!  
  a =                                                                           85!  
  
  ! Lines with varying lengths and ending in a comment after a space
  a =                                                                  75 !
  a =                                                                   76 !
  a =                                                                    77 !
  a =                                                                     78 !
  a =                                                                      79 !
  a =                                                                       80 !
  a =                                                                        81 !
  a =                                                                         82 !
  a =                                                                          83 !
  a =                                                                           84 !
  a =                                                                            85 !
  
  ! Lines with varying lengths and ending in a comment plus a space after a space
  a =                                                                 75 ! 
  a =                                                                  76 ! 
  a =                                                                   77 ! 
  a =                                                                    78 ! 
  a =                                                                     79 ! 
  a =                                                                      80 ! 
  a =                                                                       81 ! 
  a =                                                                        82 ! 
  a =                                                                         83 ! 
  a =                                                                          84 ! 
  a =                                                                           85 ! 
  
  ! Lines with varying lengths and ending in a comment plus two spaces after a space
  a =                                                                75 !  
  a =                                                                 76 !  
  a =                                                                  77 !  
  a =                                                                   78 !  
  a =                                                                    79 !  
  a =                                                                     80 !  
  a =                                                                      81 !  
  a =                                                                       82 !  
  a =                                                                        83 !  
  a =                                                                         84 !  
  a =                                                                          85 !  
  
  ! Lines with varying lengths and ending in a comment after two spaces
  a =                                                                 75  !
  a =                                                                  76  !
  a =                                                                   77  !
  a =                                                                    78  !
  a =                                                                     79  !
  a =                                                                      80  !
  a =                                                                       81  !
  a =                                                                        82  !
  a =                                                                         83  !
  a =                                                                          84  !
  a =                                                                           85  !
  
  ! Lines with varying lengths and ending in a comment plus a space after two spaces
  a =                                                                75  ! 
  a =                                                                 76  ! 
  a =                                                                  77  ! 
  a =                                                                   78  ! 
  a =                                                                    79  ! 
  a =                                                                     80  ! 
  a =                                                                      81  ! 
  a =                                                                       82  ! 
  a =                                                                        83  ! 
  a =                                                                         84  ! 
  a =                                                                          85  ! 
  
  ! Lines with varying lengths and ending in a comment plus two spaces after two spaces
  a =                                                               75  !  
  a =                                                                76  !  
  a =                                                                 77  !  
  a =                                                                  78  !  
  a =                                                                   79  !  
  a =                                                                    80  !  
  a =                                                                     81  !  
  a =                                                                      82  !  
  a =                                                                       83  !  
  a =                                                                        84  !  
  a =                                                                         85  !  
  
  ! Broken lines with varying length
  a =                                                                  75 &
       & + 75                                                             &
       + 75
  a =                                                                   76 &
       & + 76                                                              &
       + 76
  a =                                                                    77 &
       & + 77                                                               &
       + 77
  a =                                                                     78 &
       & + 78                                                                &
       + 78
  a =                                                                      79 &
       & + 79                                                                 &
       + 79
  a =                                                                       80 &
       & + 80                                                                  &
       + 80
  a =                                                                        81 &
       & + 81                                                                   &
       + 81
  a =                                                                         82 &
       & + 82                                                                    &
       + 82
  a =                                                                          83 &
       & + 83                                                                     &
       + 83
  a =                                                                           84 &
       & + 84                                                                      &
       + 84
  a =                                                                            85 &
       & + 85                                                                       &
       + 85
  
  ! Broken lines with varying length with a trailing a space
  a =                                                                 75 & 
       & + 75                                                            & 
       + 75
  a =                                                                  76 & 
       & + 76                                                             & 
       + 76
  a =                                                                   77 & 
       & + 77                                                              & 
       + 77
  a =                                                                    78 & 
       & + 78                                                               & 
       + 78
  a =                                                                     79 & 
       & + 79                                                                & 
       + 79
  a =                                                                      80 & 
       & + 80                                                                 & 
       + 80
  a =                                                                       81 & 
       & + 81                                                                  & 
       + 81
  a =                                                                        82 & 
       & + 82                                                                   & 
       + 82
  a =                                                                         83 & 
       & + 83                                                                    & 
       + 83
  a =                                                                          84 & 
       & + 84                                                                     & 
       + 84
  a =                                                                           85 & 
       & + 85                                                                      & 
       + 85
  
  ! Broken lines with varying length with two trailing a spaces
  a =                                                                75 &  
       & + 75                                                           &  
       + 75
  a =                                                                 76 &  
       & + 76                                                            &  
       + 76
  a =                                                                  77 &  
       & + 77                                                             &  
       + 77
  a =                                                                   78 &  
       & + 78                                                              &  
       + 78
  a =                                                                    79 &  
       & + 79                                                               &  
       + 79
  a =                                                                     80 &  
       & + 80                                                                &  
       + 80
  a =                                                                      81 &  
       & + 81                                                                 &  
       + 81
  a =                                                                       82 &  
       & + 82                                                                  &  
       + 82
  a =                                                                        83 &  
       & + 83                                                                   &  
       + 83
  a =                                                                         84 &  
       & + 84                                                                    &  
       + 84
  a =                                                                          85 &  
       & + 85                                                                     &  
       + 85
  
!!$  contains
!!$    
!!$    subroutine dummy
!!$      implicit none
!!$    end subroutine dummy
!!$    
end subroutine TestFreeF90
