 /*@@
   @file      Whisky_RegisterVars.F90
   @date      Thu Mar  7 20:00:56 2002
   @author    
   @desc 
   Routine to register the variables with the MoL thorn.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Whisky_Register(CCTK_ARGUMENTS)

  implicit none  

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: index1, index2, ierr

!!$  CCTK_INT :: MoLRegisterEvolvedGroup, &
!!$       MoLRegisterConstrainedGroup, MoLRegisterSaveAndRestoreGroup, &
!!$       MoLRegisterEvolved, MoLRegisterConstrained

!!$  If we''re evolving register as normal...

  ierr = 0

  call CCTK_IsFunctionAliased(ierr, "MoLRegisterEvolvedGroup")
  if (ierr == 0) then
    call CCTK_WARN(0, "The MoLRegisterEvolvedGroup function hasn't been aliased!")
  end if
  call CCTK_IsFunctionAliased(ierr, "MoLRegisterConstrainedGroup")
  if (ierr == 0) then
    call CCTK_WARN(0, "The MoLRegisterConstrainedGroup function hasn't been aliased!")
  end if
  call CCTK_IsFunctionAliased(ierr, "MoLRegisterSaveAndRestoreGroup")
  if (ierr == 0) then
    call CCTK_WARN(0, "The MoLRegisterSaveAndRestoreGroup function hasn't been aliased!")
  end if
  call CCTK_IsFunctionAliased(ierr, "MoLRegisterEvolved")
  if (ierr == 0) then
    call CCTK_WARN(0, "The MoLRegisterEvolved function hasn't been aliased!")
  end if
  call CCTK_IsFunctionAliased(ierr, "MoLRegisterConstrained")
  if (ierr == 0) then
    call CCTK_WARN(0, "The MoLRegisterConstrained function hasn't been aliased!")
  end if

  ierr = 0
  
  if (whisky_evolve .eq. 1) then

    if (CCTK_EQUALS(whisky_eos_type,"General")) then
      
      call CCTK_VarIndex(index1, "whisky::dens")
      call CCTK_VarIndex(index2, "whisky::densrhs")
      ierr = MoLRegisterEvolved(index1,index2)
      call CCTK_VarIndex(index1, "whisky::sx")
      call CCTK_VarIndex(index2, "whisky::sxrhs")
      ierr = ierr + MoLRegisterEvolved(index1,index2)
      call CCTK_VarIndex(index1, "whisky::sy")
      call CCTK_VarIndex(index2, "whisky::syrhs")
      ierr = ierr + MoLRegisterEvolved(index1,index2)
      call CCTK_VarIndex(index1, "whisky::sz")
      call CCTK_VarIndex(index2, "whisky::szrhs")
      ierr = ierr + MoLRegisterEvolved(index1,index2)
      call CCTK_VarIndex(index1, "whisky::tau")
      call CCTK_VarIndex(index2, "whisky::taurhs")
      ierr = ierr + MoLRegisterEvolved(index1,index2)
      if ( whisky_mhd_handle.eq.2 ) then
         call CCTK_VarIndex(index1, "whisky::Bnx")
         call CCTK_VarIndex(index2, "whisky::Bnxrhs")
         ierr = ierr + MoLRegisterEvolved(index1,index2)
         call CCTK_VarIndex(index1, "whisky::Bny")
         call CCTK_VarIndex(index2, "whisky::Bnyrhs")
         ierr = ierr + MoLRegisterEvolved(index1,index2)
         call CCTK_VarIndex(index1, "whisky::Bnz")
         call CCTK_VarIndex(index2, "whisky::Bnzrhs")
         ierr = ierr + MoLRegisterEvolved(index1,index2)
         if ( clean_divergence.ne.0 ) then
            call CCTK_VarIndex(index1, "whisky::divclean_psi")
            call CCTK_VarIndex(index2, "whisky::divclean_psirhs")
            ierr = ierr + MoLRegisterEvolved(index1,index2)
         end if
      end if
      if ( whisky_mhd_handle.gt.2 ) then
         call CCTK_GroupIndex(index1, "HydroBase::Avec")
         call CCTK_GroupIndex(index2, "whisky::Avecrhs")
         ierr = ierr + MoLRegisterEvolvedGroup(index1,index2)
         if ( whisky_mhd_handle.gt.3 ) then
            call CCTK_VarIndex(index1, "HydroBase::Aphi")
            call CCTK_VarIndex(index2, "whisky::Aphirhs")
            ierr = ierr + MoLRegisterEvolved(index1,index2)
         end if
      end if
 
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering conserved variables with MoL")
      end if
      
!!$      call CCTK_GroupIndex(index1, "whisky::whisky_cons_var")
!!$      call CCTK_GroupIndex(index2, "whisky::whisky_rhs_var")
!!$      ierr = MoLRegisterEvolvedGroup(index1,index2)
!!$      
!!$      if (ierr .ne. 0) then
!!$
!!$        call CCTK_WARN(0, "Error registering conserved variables with MoL")
!!$
!!$      end if
      
    else if (CCTK_EQUALS(whisky_eos_type,"Polytype")) then
      
      call CCTK_VarIndex(index1, "whisky::dens")
      call CCTK_VarIndex(index2, "whisky::densrhs")
      ierr = MoLRegisterEvolved(index1, index2)
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering conserved D with MoL")
      end if
      
      call CCTK_VarIndex(index1, "whisky::sx")
      call CCTK_VarIndex(index2, "whisky::sxrhs")
      ierr = MoLRegisterEvolved(index1, index2)
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering conserved Sx with MoL")
      end if

      call CCTK_VarIndex(index1, "whisky::sy")
      call CCTK_VarIndex(index2, "whisky::syrhs")
      ierr = MoLRegisterEvolved(index1, index2)
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering conserved Sy with MoL")
      end if
      
      call CCTK_VarIndex(index1, "whisky::sz")
      call CCTK_VarIndex(index2, "whisky::szrhs")
      ierr = MoLRegisterEvolved(index1, index2)
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering conserved Sz with MoL")
      end if
      
      call CCTK_VarIndex(index1, "whisky::tau")
      ierr = MoLRegisterConstrained(index1)
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering conserved tau with MoL")
      end if

      if ( whisky_mhd_handle.eq.2 ) then
         call CCTK_VarIndex(index1, "whisky::Bnx")
         call CCTK_VarIndex(index2, "whisky::Bnxrhs")
         ierr = MoLRegisterEvolved(index1,index2)
         if (ierr .ne. 0) then
           call CCTK_WARN(0, "Error registering conserved Bnx with MoL")
         end if

         call CCTK_VarIndex(index1, "whisky::Bny")
         call CCTK_VarIndex(index2, "whisky::Bnyrhs")
         ierr = MoLRegisterEvolved(index1,index2)
         if (ierr .ne. 0) then
           call CCTK_WARN(0, "Error registering conserved Bny with MoL")
         end if

         call CCTK_VarIndex(index1, "whisky::Bnz")
         call CCTK_VarIndex(index2, "whisky::Bnzrhs")
         ierr = MoLRegisterEvolved(index1,index2)
         if (ierr .ne. 0) then
           call CCTK_WARN(0, "Error registering conserved Bnz with MoL")
         end if

         if ( clean_divergence.ne.0 ) then
            call CCTK_VarIndex(index1, "whisky::divclean_psi") 
            call CCTK_VarIndex(index2, "whisky::divclean_psirhs") 
            ierr = MoLRegisterEvolved(index1,index2)
            if ( ierr.ne.0 ) then
               call CCTK_WARN(0, "Error registering divclean_psi with MoL")
            end if
         end if 
      end if
      if ( whisky_mhd_handle.gt.2 ) then
         call CCTK_GroupIndex(index1, "HydroBase::Avec")
         call CCTK_GroupIndex(index2, "whisky::Avecrhs")
         ierr = MoLRegisterEvolvedGroup(index1,index2)
         if ( ierr.ne.0 ) then
            call CCTK_WARN(0, "Error registering Avec with MoL")
         end if
         if ( whisky_mhd_handle.gt.3 ) then
            call CCTK_VarIndex(index1, "HydroBase::Aphi")
            call CCTK_VarIndex(index2, "whisky::Aphirhs")
            ierr = MoLRegisterEvolved(index1,index2)
            if ( ierr.ne.0 ) then
               call CCTK_WARN(0, "Error registering Aphi with MoL")
            end if
         end if
      end if
      
    else
      call CCTK_WARN(0, "Don't recognize the type of EOS!")
    end if
    

!!$    call CCTK_VarIndex(index1, "whisky::dens")
!!$    call CCTK_VarIndex(index2, "whisky::densrhs")
!!$    call MoL_RegisterVar(ierr, index1,index2)
!!$    call CCTK_VarIndex(index1, "whisky::sx")
!!$    call CCTK_VarIndex(index2, "whisky::sxrhs")
!!$    call MoL_RegisterVar(ierr, index1,index2)
!!$    call CCTK_VarIndex(index1, "whisky::sy")
!!$    call CCTK_VarIndex(index2, "whisky::syrhs")
!!$    call MoL_RegisterVar(ierr, index1,index2)
!!$    call CCTK_VarIndex(index1, "whisky::sz")
!!$    call CCTK_VarIndex(index2, "whisky::szrhs")
!!$    call MoL_RegisterVar(ierr, index1,index2)
!!$    
!!$    if (CCTK_EQUALS(whisky_eos_type,"General")) then
!!$      
!!$      call CCTK_VarIndex(index1, "whisky::tau")
!!$      call CCTK_VarIndex(index2, "whisky::taurhs")
!!$      call MoL_RegisterVar(ierr, index1,index2)
!!$      
!!$    else if (CCTK_EQUALS(whisky_eos_type,"Polytype")) then
!!$      
!!$      call CCTK_VarIndex(index1, "whisky::tau")
!!$      call MoL_RegisterPrimitive(ierr, index1)
!!$      
!!$    else
!!$      
!!$      call CCTK_WARN(0, "Don't recognize the type of EOS!")
!!$      
!!$    end if
!!$    

    call CCTK_VarIndex(index1, "HydroBase::rho")
    ierr = MoLRegisterConstrained(index1)
    call CCTK_VarIndex(index1, "HydroBase::press")
    ierr = ierr + MoLRegisterConstrained(index1)
    call CCTK_VarIndex(index1, "HydroBase::eps")
    ierr = ierr + MoLRegisterConstrained(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering primitive variables with MoL")
    end if

    call CCTK_GroupIndex(index1, "Whisky::whisky_aux_scalar_var")
    ierr = MoLRegisterConstrainedGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering auxilliary variables with MoL")
    end if

    call CCTK_GroupIndex(index1, "HydroBase::vel")
    ierr = MoLRegisterConstrainedGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering primitive velocity variable with MoL")
    end if

    if ( whisky_mhd_handle.gt.0 ) then

      call CCTK_GroupIndex(index1, "HydroBase::Bvec")
      ierr = MoLRegisterConstrainedGroup(index1)
      if (ierr .ne. 0) then
         call CCTK_WARN(0, "Error registering primitive magnetic variables with MoL")
      end if

      if ( monitor_divB.ne.0 ) then
          call CCTK_GroupIndex(index1, "Whisky::divB")
          ierr = MoLRegisterConstrainedGroup(index1)
          if (ierr.ne.0) then
             call CCTK_WARN(0, "Error registering divB with MoL")
          end if
      end if
      if ( calculate_b2_comoving.ne.0 ) then
          call CCTK_GroupIndex(index1, "Whisky::b2_comoving")
          ierr = MoLRegisterConstrainedGroup(index1)
          if (ierr.ne.0) then
             call CCTK_WARN(0, "Error registering b2_comoving with MoL")
          end if
      end if

      if ( whisky_mhd_handle.gt.2 ) then
         call CCTK_GroupIndex(index1, "Whisky::whisky_Bfield_var")
         ierr = MoLRegisterConstrainedGroup(index1)
         if (ierr.ne.0) then
            call CCTK_WARN(0, "Error registering Bn[xyz] as constrained with MoL" )
         end if
      end if

    end if

!!$    call CCTK_VarIndex(index1, "whisky::rho")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::velx")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::vely")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::velz")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::press")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::eps")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::w_lorentz")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    

    call CCTK_GroupIndex(index1, "admbase::lapse")
    ierr = MoLRegisterSaveAndRestoreGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering lapse with MoL")
    end if

    call CCTK_GroupIndex(index1, "admbase::metric")
    ierr = MoLRegisterSaveAndRestoreGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering metric with MoL")
    end if

    call CCTK_GroupIndex(index1, "admbase::curv")
    ierr = MoLRegisterSaveAndRestoreGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering curvature with MoL")
    end if

!!$    call CCTK_VarIndex(index1, "ADMBase::alp")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::gxx")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::gxy")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::gxz")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::gyy")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::gyz")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::gzz")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::kxx")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::kxy")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::kxz")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::kyy")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::kyz")
!!$    call MoL_RegisterDepends(ierr, index1)
!!$    call CCTK_VarIndex(index1, "ADMBase::kzz")
!!$    call MoL_RegisterDepends(ierr, index1)
    
    if (.not.CCTK_EQUALS(initial_shift,"none")) then

      if (CCTK_EQUALS(shift_evolution_method,"Comoving")) then

        call CCTK_GroupIndex(index1, "admbase::shift")
        ierr = MoLRegisterConstrainedGroup(index1)
        if (ierr .ne. 0) then
          call CCTK_WARN(0, "Error registering shift with MoL")
        end if

        call CCTK_GroupIndex(index1, "Whisky::whisky_coords")
        call CCTK_GroupIndex(index2, "Whisky::whisky_coords_rhs")
        ierr = MoLRegisterEvolvedGroup(index1, index2)
        if (ierr .ne. 0) then
          call CCTK_WARN(0, "Error registering comoving coords with MoL")
        end if

      else
      
        call CCTK_GroupIndex(index1, "admbase::shift")
        ierr = MoLRegisterSaveAndRestoreGroup(index1)
        if (ierr .ne. 0) then
          call CCTK_WARN(0, "Error registering shift with MoL")
        end if
      
      end if
      
!!$      call CCTK_VarIndex(index1, "ADMBase::betax")
!!$      call MoL_RegisterDepends(ierr, index1)
!!$      call CCTK_VarIndex(index1, "ADMBase::betay")
!!$      call MoL_RegisterDepends(ierr, index1)
!!$      call CCTK_VarIndex(index1, "ADMBase::betaz")
!!$      call MoL_RegisterDepends(ierr, index1)
      
    end if
    
    if (evolve_tracer .ne. 0) then

      call CCTK_GroupIndex(index1, "whisky::whisky_cons_tracers")
      call CCTK_GroupIndex(index2, "whisky::whisky_tracer_rhs")
      ierr = MoLRegisterEvolvedGroup(index1, index2)
      if (ierr .ne. 0) then
        call CCTK_WARN(0, "Error registering tracers with MoL")
      end if
      
    end if
   
    if (number_of_particles .gt. 0) then

      call CCTK_GroupIndex(index1, "whisky::particles")
      call CCTK_GroupIndex(index2, "whisky::particle_rhs")
      ierr = MoLRegisterEvolvedGroup(index1,index2)

    end if
    
  else !! if (whisky_evolve .ne. 1), just register as constrained (why?)

    call CCTK_GroupIndex(index1, "whisky::whisky_cons_scalar_var")
    ierr = MoLRegisterConstrainedGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering conserved variables with MoL")
    end if

    call CCTK_GroupIndex(index1, "whisky::whisky_cons_vector_var")
    ierr = MoLRegisterConstrainedGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering conserved variables with MoL")
    end if

    if ( whisky_mhd_handle.ge.2 ) then
       call CCTK_GroupIndex(index1, "whisky::whisky_Bfield_var")
       ierr = MoLRegisterConstrainedGroup(index1)
       if (ierr .ne. 0) then
         call CCTK_WARN(0, "Error registering conserved magnetic field variables with MoL")
       end if
    end if
    if ( whisky_mhd_handle.ge.3 ) then
       call CCTK_GroupIndex(index1, "HydroBase::Avec")
       ierr = MoLRegisterConstrainedGroup(index1)
       if (ierr .ne. 0) then
         call CCTK_WARN(0, "Error registering Avec with MoL")
       end if
    end if
    if ( whisky_mhd_handle.ge.4 ) then
       call CCTK_GroupIndex(index1, "HydroBase::Aphi")
       ierr = MoLRegisterConstrainedGroup(index1)
       if (ierr .ne. 0) then
         call CCTK_WARN(0, "Error registering Aphi with MoL")
       end if
    end if


    call CCTK_VarIndex(index1, "HydroBase::rho")
    ierr = MoLRegisterConstrained(index1)
    call CCTK_VarIndex(index1, "HydroBase::press")
    ierr = MoLRegisterConstrained(index1)
    call CCTK_VarIndex(index1, "HydroBase::eps")
    ierr = MoLRegisterConstrained(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering primitive variables with MoL")
    end if

    call CCTK_GroupIndex(index1, "Whisky::whisky_aux_scalar_var")
    ierr = MoLRegisterConstrainedGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering auxilliary variables with MoL")
    end if

    call CCTK_GroupIndex(index1, "HydroBase::vel")
    ierr = MoLRegisterConstrainedGroup(index1)
    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Error registering primitive variables with MoL")
    end if

    if ( whisky_mhd_handle.gt.0 ) then
       call CCTK_GroupIndex(index1, "HydroBase::Bvec")
       ierr = MoLRegisterConstrainedGroup(index1)
       if (ierr .ne. 0) then
         call CCTK_WARN(0, "Error registering primitive magnetic variables with MoL")
       end if
    end if
    

!!$    call CCTK_VarIndex(index1, "whisky::dens")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::sx")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::sy")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::sz")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::tau")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    
!!$    call CCTK_VarIndex(index1, "whisky::rho")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::velx")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::vely")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::velz")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::press")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::eps")
!!$    call MoL_RegisterPrimitive(ierr, index1)
!!$    call CCTK_VarIndex(index1, "whisky::w_lorentz")
!!$    call MoL_RegisterPrimitive(ierr, index1)
    
  end if

end subroutine Whisky_Register
