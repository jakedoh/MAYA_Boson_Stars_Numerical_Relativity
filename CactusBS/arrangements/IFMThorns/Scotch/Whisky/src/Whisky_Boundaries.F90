 /*@@
   @file      Whisky_Boundaries.F90
   @date      Sat Jan 26 01:01:14 2002
   @author    
   @desc 
   The two routines for dealing with boundary conditions.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "util_Table.h"
#include "Whisky_Utils.h"

 /*@@
   @routine    Whisky_InitSymBound
   @date       Sat Jan 26 01:03:04 2002
   @author     Ian Hawke
   @desc 
   Sets up the symmetries at the boundaries of the hydrodynamical variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Direct translation of routines from GR3D, GRAstro_Hydro, 
   written by Mark Miller, or WaveToy routines, or...
   @endhistory 

@@*/

subroutine Whisky_InitSymBound(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer, dimension(3) :: sym
  integer :: ierr
  integer :: itracer

  character(len=100) tracername
  character(len=100) tracerindex

  sym(1) = 1
  sym(2) = 1
  sym(3) = 1
  
  if (sync_primitives .ne. 0) then
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::rho")
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::press")
    call SetCartSymVN(ierr, cctkGH, sym, "Whisky::w_lorentz")
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::eps")
  end if
  call SetCartSymVN(ierr, cctkGH, sym, "Whisky::dens")
  call SetCartSymVN(ierr, cctkGH, sym, "Whisky::tau")
  if (sync_primitives .eq. 0) then ! only apply BC if primitives are NOT synced
    call SetCartSymVN(ierr, cctkGH, sym, "spacemask::space_mask");
  end if

!!$ handle multiple tracer variables
  if(evolve_tracer.eq.1) then
     call SetCartSymGN(ierr, cctkGH, sym, "Whisky::whisky_tracers")
     call SetCartSymGN(ierr, cctkGH, sym, "Whisky::whisky_cons_tracers")
  endif 
  
  sym(1) = -1
  sym(2) = 1
  sym(3) = 1
  
  if (sync_primitives .ne. 0) then 
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::velx")
  end if
  call SetCartSymVN(ierr, cctkGH, sym, "Whisky::sx")
  
  sym(1) = 1
  sym(2) = -1
  sym(3) = 1
  
  if (sync_primitives .ne. 0) then 
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::vely")
  end if
  call SetCartSymVN(ierr, cctkGH, sym, "Whisky::sy")
  
  sym(1) = 1
  sym(2) = 1
  sym(3) = -1
  
  if (sync_primitives .ne. 0) then 
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::velz")
  end if
  call SetCartSymVN(ierr, cctkGH, sym, "Whisky::sz")

  if(whisky_mhd_handle.gt.0) then

    sym(1) = -1
    sym(2) = 1
    sym(3) = 1
    if (whisky_mhd_handle.gt.1) then
       call SetCartSymVN(ierr, cctkGH, sym, "Whisky::Bnx")
    end if
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Bvec[0]")

    sym(1) = 1
    sym(2) = -1
    sym(3) = 1
    if (whisky_mhd_handle.gt.1) then
       call SetCartSymVN(ierr, cctkGH, sym, "Whisky::Bny")
    end if
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Bvec[1]")

    sym(1) = 1
    sym(2) = 1
    sym(3) = -1
    if (whisky_mhd_handle.gt.1) then
       call SetCartSymVN(ierr, cctkGH, sym, "Whisky::Bnz")
    end if
    call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Bvec[2]")

    if (clean_divergence.ne.0) then

       sym(1) = 1
       sym(2) = 1
       sym(3) = 1
       call SetCartSymVN(ierr, cctkGH, sym, "Whisky::divclean_psi")

    end if

    if (monitor_divB.ne.0) then

       sym(1) = 1
       sym(2) = 1
       sym(3) = 1
       call SetCartSymVN(ierr, cctkGH, sym, "Whisky::divB")

    end if

    if (whisky_mhd_handle.gt.2) then

       sym(1) = -1
       sym(2) = 1
       sym(3) = 1
       call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Avec[0]")

       sym(1) = 1
       sym(2) = -1 
       sym(3) = 1
       call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Avec[1]")

       sym(1) = 1
       sym(2) = 1
       sym(3) = -1
       call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Avec[2]")

       if (whisky_mhd_handle.gt.3) then

          sym(1) = 1
          sym(2) = 1
          sym(3) = 1
          call SetCartSymVN(ierr, cctkGH, sym, "HydroBase::Aphi")

       end if

    end if

  end if
  
end subroutine Whisky_InitSymBound


 /*@@
   @routine    Whisky_PrimitiveBoundaries
   @date       Fri Sep 11 11:41:55 EDT 2009
   @author     
   @desc 
   Calls the appropriate boundary routines for the primitives
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_SelectPrimitiveBoundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: ierr = -1 ! make sure we throw an error if an unknown boundary condition is chosen

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  CCTK_INT, parameter :: izero=0
  
!!$Flat boundaries if required  

  if (CCTK_EQUALS(bound,"flat")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_aux_scalar_var", "Flat")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::rho", "Flat")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::press", "Flat")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::eps", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::vel", "Flat")
    if (whisky_mhd_handle.gt.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Bvec", "Flat")
    end if
    if (whisky_mhd_handle.gt.2) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Avec", "Flat")
    end if
    if (whisky_mhd_handle.gt.3) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Aphi", "Flat")
    end if
  endif

  if (CCTK_EQUALS(bound,"none")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_aux_scalar_var", "None")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::rho", "None")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::press", "None")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::eps", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::vel", "None")
    if (whisky_mhd_handle.gt.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Bvec", "None")
    end if
    if (whisky_mhd_handle.gt.2) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Avec", "None")
    end if
    if (whisky_mhd_handle.gt.3) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Aphi", "None")
    end if
  end if

  if (CCTK_EQUALS(bound,"static")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_aux_scalar_var", "Static")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::rho", "Static")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::press", "Static")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::eps", "Static")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::vel", "Static")
    if (whisky_mhd_handle.gt.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Bvec", "Static")
    end if
    if (whisky_mhd_handle.gt.2) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Avec", "Static")
    end if
    if (whisky_mhd_handle.gt.3) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Aphi", "Static")
    end if
  end if
  
  
  !!---------------------HACK
  if (CCTK_EQUALS(bound,"hybrid")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_aux_scalar_var", "custom")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::rho", "custom")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::press", "custom")
    ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::eps", "custom")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "HydroBase::vel", "custom")
    if (whisky_mhd_handle.gt.0) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Bvec", "custom")
    end if
    if (whisky_mhd_handle.gt.2) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Avec", "custom")
    end if
    if (whisky_mhd_handle.gt.3) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "HydroBase::Aphi", "custom")
    end if
  end if
  
  
  
  !!---------------------End Hack---

  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
    return ! NOTREACHED
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine Whisky_SelectPrimitiveBoundaries

 /*@@
   @routine    Whisky_Boundaries
   @date       Sat Jan 26 01:04:04 2002
   @author     
   @desc 
   Calls the appropriate boundary routines
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_Boundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3) :: sw
  integer :: ierr = -1 ! make sure we throw an error if an unknown boundary condition is chosen
  CCTK_REAL :: pi = 3.141569d0
  integer :: i,j,k

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  CCTK_INT, parameter :: izero=0
  
  sw = Whisky_stencil

!!$Symmetry boundaries first  

!!$  call CartSymGN(ierr,cctkGH,"Whisky::Whisky_cons_var")
!!$  call CartSymVN(ierr,cctkGH,"HydroBase::rho")
!!$  call CartSymVN(ierr,cctkGH,"HydroBase::press")
!!$  call CartSymVN(ierr,cctkGH,"HydroBase::eps")
!!$  call CartSymGN(ierr,cctkGH,"HydroBase::vel")

  
!!$Flat boundaries if required  

  if (CCTK_EQUALS(bound,"flat")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_scalar_var", "Flat")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_vector_var", "Flat")
    if (whisky_mhd_handle.gt.1) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::whisky_bfield_var", "Flat")
       if (clean_divergence.ne.0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
             "Whisky::divclean_psi", "Flat")
       end if
       if (monitor_divB .ne. 0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::divB", "Flat")
       end if
       if (whisky_mhd_handle.gt.2) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Avec", "Flat")
       end if
       if (whisky_mhd_handle.gt.3) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Aphi", "Flat")
       end if
    end if


    if (sync_primitives .ne. 0) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "Whisky::Whisky_aux_scalar_var", "Flat")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::rho", "Flat")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::press", "Flat")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::eps", "Flat")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::vel", "Flat")
      if (whisky_mhd_handle.gt.0) then
         ierr = ierr +Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
              "HydroBase::Bvec", "Flat")
      end if
    end if

    if(evolve_tracer .ne. 0) then 
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_tracers", "Flat")
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_cons_tracers", "Flat")
    endif


!!$    call BndFlatGN(ierr,cctkGH,sw,"Whisky::Whisky_cons_var")
!!$    call BndFlatVN(ierr,cctkGH,sw,"HydroBase::rho")
!!$    call BndFlatVN(ierr,cctkGH,sw,"HydroBase::press")
!!$    call BndFlatVN(ierr,cctkGH,sw,"HydroBase::eps")
!!$    call BndFlatGN(ierr,cctkGH,sw,"HydroBase::vel")
  endif

  if (CCTK_EQUALS(bound,"none")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_scalar_var", "None")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_vector_var", "None")
    if (whisky_mhd_handle.gt.1) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::whisky_bfield_var", "None")
       if (clean_divergence.ne.0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
             "Whisky::divclean_psi", "None")
       end if
       if (monitor_divB .ne. 0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::divB", "None")
       end if
       if (whisky_mhd_handle.gt.2) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Avec", "None")
       end if
       if (whisky_mhd_handle.gt.3) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Aphi", "None")
       end if
    end if

    if (sync_primitives .ne. 0) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "Whisky::Whisky_aux_scalar_var", "None")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::rho", "None")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::press", "None")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::eps", "None")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::vel", "None")
      if (whisky_mhd_handle.gt.0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
              "HydroBase::Bvec", "None")
      end if
    end if
!!$    call  BndStaticGN(ierr,cctkGH,sw,"Whisky::Whisky_cons_var")
!!$    call  BndStaticGN(ierr,cctkGH,sw,"Whisky::Whisky_prim_var")

    if(evolve_tracer .ne. 0) then 
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_tracers", "None")
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_cons_tracers", "None")
    end if

  end if
  ! hack
  if (CCTK_EQUALS(bound,"hybrid")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_scalar_var", "custom")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_vector_var", "custom")
    if (whisky_mhd_handle.gt.1) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::whisky_bfield_var", "custom")
       if (clean_divergence.ne.0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
             "Whisky::divclean_psi", "custom")
       end if
       if (monitor_divB .ne. 0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::divB", "custom")
       end if
       if (whisky_mhd_handle.gt.2) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Avec", "custom")
       end if
       if (whisky_mhd_handle.gt.3) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Aphi", "custom")
       end if
    end if


    if (sync_primitives .ne. 0) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "Whisky::Whisky_aux_scalar_var", "custom")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::rho", "custom")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::press", "custom")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::eps", "custom")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::vel", "custom")
      if (whisky_mhd_handle.gt.0) then
         ierr = ierr +Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
              "HydroBase::Bvec", "custom")
      end if
    end if

    if(evolve_tracer .ne. 0) then 
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_tracers", "custom")
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_cons_tracers", "custom")
    endif

  endif
  
  !----------------------
  !! Hack-------------
  if (CCTK_EQUALS(bound,"static")) then
    ierr = izero
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_scalar_var", "Static")
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
         "Whisky::Whisky_cons_vector_var", "Static")
    if (whisky_mhd_handle.gt.1) then
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::whisky_bfield_var", "Static")
       if (clean_divergence.ne.0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
             "Whisky::divclean_psi", "Static")
       end if
       if (monitor_divB .ne. 0) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::divB", "Static")
       end if
       if (whisky_mhd_handle.gt.2) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Avec", "Static")
       end if
       if (whisky_mhd_handle.gt.3) then
          ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
               "HydroBase::Aphi", "Static")
       end if
    end if

    if (sync_primitives .ne. 0) then
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "Whisky::Whisky_aux_scalar_var", "Static")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::rho", "Static")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::press", "Static")
      ierr = ierr + Boundary_SelectVarForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::eps", "Static")
      ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
           "HydroBase::vel", "Static")
      if (whisky_mhd_handle.gt.0) then
         ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
              "HydroBase::Bvec", "Static")
      end if
    end if

    if(evolve_tracer .ne. 0) then 
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_tracers", "Static")
       ierr = ierr + Boundary_SelectGroupForBC(cctkGH, faces, whisky_stencil, -ione, &
            "Whisky::Whisky_cons_tracers", "Static")
    endif

!!$    call  BndStaticGN(ierr,cctkGH,sw,"Whisky::Whisky_cons_var")
!!$    call  BndStaticGN(ierr,cctkGH,sw,"Whisky::Whisky_prim_var")
  end if

  if (sync_primitives .eq. 0) then ! only apply BC if primitives are NOT synced
    ierr = ierr + Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, whisky_stencil, -ione, &
                      "SpaceMask::space_mask_group", "none");
  end if

  if (CCTK_EQUALS(bound,"scalar")) then
    call CCTK_WARN(0, "Until somebody uses this I see no reason to support it")
!!$    call  BndScalarGN(ierr,cctkGH,sw,zero,"Whisky::Whisky_cons_var")
!!$    call  BndScalarGN(ierr,cctkGH,sw,zero,"Whisky::Whisky_prim_var")
  end if

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine Whisky_Boundaries


 /*@@
   @routine    Whisky_OutflowBoundaries
   @date       Tue May 17 16:58:06 2005
   @author     Ian Hawke
   @desc 
   Set outflow boundaries over only part of the domain.
   This is designed to be used with GZPatchSystem.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_OutflowBoundaries(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer, dimension(3) :: sw
  integer :: i,j,k
  
  CCTK_REAL, dimension(3) :: posn

  CCTK_REAL :: det

  sw = Whisky_stencil

  if (r(1,1,1) < r(1,1,cctk_lsh(3))) then

    if (cctk_bbox(6) .ne. 0) then

      do k = cctk_lsh(3) - sw(3), cctk_lsh(3)
        !$OMP PARALLEL DO PRIVATE (i,j, posn, det)
        do j = 1, cctk_lsh(2)
          do i = 1, cctk_lsh(1)
            
            posn(1) = x(i,j,k)
            posn(2) = y(i,j,k)
            posn(3) = z(i,j,k)
            
            if (dot_product(outflowboundary_normal, posn) > 0.d0) then
              
              rho(i,j,k)        = rho(i,j,k-1)
              velx(i,j,k)       = velx(i,j,k-1)
              vely(i,j,k)       = vely(i,j,k-1)
              velz(i,j,k)       = velz(i,j,k-1)
              eps(i,j,k)        = eps(i,j,k-1)
              press(i,j,k)      = press(i,j,k-1)
              w_lorentz(i,j,k)  = w_lorentz(i,j,k-1)
              if ( whisky_mhd_handle.gt.0 ) then
                 call CCTK_WARN(1,'We are putting junk into the magnetic field at the boundaries in Whisky_OutflowBoundaries.')
                 bvec(i,j,k,1) = bvec(i,j,k-1,1)
                 bvec(i,j,k,2) = bvec(i,j,k-1,2)
                 bvec(i,j,k,3) = bvec(i,j,k-1,3)
              end if
              if ( whisky_mhd_handle.gt.2 ) then
                 call CCTK_WARN(1, 'We are putting junk into Avec at the outer boundaries in Whisky_OutflowBoundaries.')
                 Avec(i,j,k,1) = Avec(i,j,k-1,1)
                 Avec(i,j,k,2) = Avec(i,j,k-1,2)
                 Avec(i,j,k,3) = Avec(i,j,k-1,3)
              end if
              if ( whisky_mhd_handle.gt.3 ) then
                 call CCTK_WARN(1, 'We are putting junk into Avec at the outer boundaries in Whisky_OutflowBoundaries.')
                 Aphi(i,j,k) = Aphi(i,j,k-1)
              end if
              
              det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if ( whisky_mhd_handle.gt.1 ) then
                 call prim2conmhd(whisky_eos_handle,gxx(i,j,k),&
                      gxy(i,j,k),gxz(i,j,k),&
                      gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                      det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                      tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                      eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),bnx(i,j,k),&
                      bny(i,j,k),bnz(i,j,k),bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
              else
                 call prim2con(whisky_eos_handle,gxx(i,j,k),&
                      gxy(i,j,k),gxz(i,j,k),&
                      gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                      det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                      tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                      eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
              end if
 
              
            end if
          
          end do
        end do
        !$OMP END PARALLEL DO
      end do
      
    end if
  
  else

    if (cctk_bbox(5) .ne. 0) then

      do k = sw(3), 1, -1
        !$OMP PARALLEL DO PRIVATE (i,j, posn, det)
        do j = 1, cctk_lsh(2)
          do i = 1, cctk_lsh(1)
            
            posn(1) = x(i,j,k)
            posn(2) = y(i,j,k)
            posn(3) = z(i,j,k)
            
            if (dot_product(outflowboundary_normal, posn) > 0.d0) then
              
              rho(i,j,k)        = rho(i,j,k+1)
              velx(i,j,k)       = velx(i,j,k+1)
              vely(i,j,k)       = vely(i,j,k+1)
              velz(i,j,k)       = velz(i,j,k+1)
              eps(i,j,k)        = eps(i,j,k+1)
              press(i,j,k)      = press(i,j,k+1)
              w_lorentz(i,j,k)  = w_lorentz(i,j,k+1)
              if ( whisky_mhd_handle.gt.1 ) then
                 call CCTK_WARN(1,'We are putting junk into the magnetic field at the boundaries in Whisky_OutflowBoundaries.')
                 bvec(i,j,k,1) = bvec(i,j,k+1,1)
                 bvec(i,j,k,2) = bvec(i,j,k+1,2)
                 bvec(i,j,k,3) = bvec(i,j,k+1,3)
              end if
              if ( whisky_mhd_handle.gt.2 ) then
                 call CCTK_WARN(1, 'We are putting junk into Avec at the outer boundaries in Whisky_OutflowBoundaries.')
                 Avec(i,j,k,1) = Avec(i,j,k+1,1)
                 Avec(i,j,k,2) = Avec(i,j,k+1,2)
                 Avec(i,j,k,3) = Avec(i,j,k+1,3)
              end if
              if ( whisky_mhd_handle.gt.3 ) then
                 call CCTK_WARN(1, 'We are putting junk into Avec at the outer boundaries in Whisky_OutflowBoundaries.')
                 Aphi(i,j,k) = Aphi(i,j,k+1)
              end if
              
              
              det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
              if ( whisky_mhd_handle.gt.1 ) then
                 call prim2conmhd(whisky_eos_handle,gxx(i,j,k),&
                      gxy(i,j,k),gxz(i,j,k),&
                      gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                      det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                      tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                      eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),bnx(i,j,k),&
                      bny(i,j,k),bnz(i,j,k),bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
              else
                 call prim2con(whisky_eos_handle,gxx(i,j,k),&
                      gxy(i,j,k),gxz(i,j,k),&
                      gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                      det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                      tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                      eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
              end if
              
            end if
          
          end do
        end do
        !$OMP END PARALLEL DO
      end do
      
    end if

  end if
  
end subroutine Whisky_OutflowBoundaries

