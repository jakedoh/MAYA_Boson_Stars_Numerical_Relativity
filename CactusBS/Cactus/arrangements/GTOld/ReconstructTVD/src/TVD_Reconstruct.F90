 /*@@
   @file      TVD_Reconstruct.F90
   @date      Sat Jan 26 02:13:25 2002
   @author    
   @desc 
   Wrapper routine to perform the reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

 /*@@
   @routine    TVD_Reconstruct
   @date       Sat Jan 26 02:13:47 2002
   @author     Luca Baiotti, Ian Hawke
   @desc 
   A wrapper routine to do reconstruction. Currently just does
   TVD on the primitive variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine TVD_Reconstruct(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: nx, ny, nz, i, j, k, itracer

  logical, dimension(:,:,:), allocatable :: trivial_rp
!!$   logical, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz
  CCTK_INT, dimension(3) :: excision_descriptors

  CCTK_INT :: ierr

  CCTK_REAL :: local_min_tracer

  allocate(trivial_rp(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Allocation problems with trivial_rp")
  end if
  
    excision_descriptors(1)=-1
    excision_descriptors(2)=-1
    excision_descriptors(3)=-1

  call SpaceMask_GetTypeBits(type_bitsx, "Hydro_RiemannProblemX")
  call SpaceMask_GetStateBits(trivialx, "Hydro_RiemannProblemX", &
       &"trivial")
  call SpaceMask_GetStateBits(not_trivialx, "Hydro_RiemannProblemX", &
       &"not_trivial")
  call SpaceMask_GetTypeBits(type_bitsy, "Hydro_RiemannProblemY")
  call SpaceMask_GetStateBits(trivialy, "Hydro_RiemannProblemY", &
       &"trivial")
  call SpaceMask_GetStateBits(not_trivialy, "Hydro_RiemannProblemY", &
       &"not_trivial")
  call SpaceMask_GetTypeBits(type_bitsz, "Hydro_RiemannProblemZ")
  call SpaceMask_GetStateBits(trivialz, "Hydro_RiemannProblemZ", &
       &"trivial")
  call SpaceMask_GetStateBits(not_trivialz, "Hydro_RiemannProblemZ", &
       &"not_trivial")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

!!$  Currently only option is reconstruction on primitive variables.
!!$  Should change this.

!!$ Initialize variables that store reconstructed quantities
!!$ NOTE: conservative reconstruction is broken since its con2prim does not have
!!$       starting guesses for the iteration

  rhoplus = 0.0d0
  rhominus = 0.0d0
  epsplus = 0.0d0
  epsminus = 0.0d0
  velxplus = 0.0d0
  velxminus = 0.0d0
  velyplus = 0.0d0
  velyminus = 0.0d0
  velzplus = 0.0d0
  velzminus = 0.0d0
  if (whisky_mhd_handle.gt.1) then
    Bxplus = 0.0d0
    Bxminus = 0.0d0
    Byplus = 0.0d0
    Byminus = 0.0d0
    Bzplus = 0.0d0
    Bzminus = 0.0d0
    if (clean_divergence.ne.0) then
       divclean_psiplus = 0.0d0
       divclean_psiminus = 0.0d0
    end if

  end if

  if (evolve_tracer .ne. 0) then
     tracerplus = 0.0d0
     tracerminus = 0.0d0
  endif

     if (evolve_tracer .ne. 0) then
        do itracer=1,number_of_tracers
           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                tracer(:,:,:,itracer), tracerplus(:,:,:,itracer), tracerminus(:,:,:,itracer), & 
                trivial_rp, space_mask, excision_descriptors)        
        enddo
     end if

    if (CCTK_EQUALS(recon_vars,"primitive")) then

      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           rho, rhoplus, rhominus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           velx, velxplus, velxminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           vely, velyplus, velyminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           velz, velzplus, velzminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           eps, epsplus, epsminus, trivial_rp, space_mask, &
           excision_descriptors)

      if (whisky_mhd_handle.gt.1) then
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Bvec(:,:,:,1), Bxplus, Bxminus, trivial_rp, space_mask, &
             excision_descriptors)
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Bvec(:,:,:,2), Byplus, Byminus, trivial_rp, space_mask, &
             excision_descriptors)
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Bvec(:,:,:,3), Bzplus, Bzminus, trivial_rp, space_mask, &
             excision_descriptors)

        if (clean_divergence.ne.0) then

           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                divclean_psi, divclean_psiplus, divclean_psiminus, &
                trivial_rp, space_mask, excision_descriptors)

        end if

      end if

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if (trivial_rp(i,j,k)) then
              if (flux_direction == 1) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsx, trivialx)
              else if (flux_direction == 2) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsy, trivialy)
              else if (flux_direction == 3) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsz, trivialz)
              end if
            else
              if (flux_direction == 1) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsx, not_trivialx)
              else if (flux_direction == 2) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsy, not_trivialy)
              else if (flux_direction == 3) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsz, not_trivialz)
              end if
            end if
          end do
        end do
      end do
    else if (CCTK_EQUALS(recon_vars,"conservative")) then
      call CCTK_WARN(0, "Reconstruction of conservatives is broken since no \
                         initial guess for the primitive plus/minus variables \
                         is provided")
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           dens, densplus, densminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           sx, sxplus, sxminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           sy, syplus, syminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           sz, szplus, szminus, trivial_rp, space_mask, &
           excision_descriptors)
      call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
           tau, tauplus, tauminus, trivial_rp, space_mask, &
           excision_descriptors)
      if (whisky_mhd_handle.gt.1) then
        call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
             Bnx, Bnxplus, Bnxminus, trivial_rp, space_mask, &
             excision_descriptors)
        call tvdreconstruct(ny, ny, nz, xoffset, yoffset, zoffset, &
             Bny, Bnyplus, Bnyminus, trivial_rp, space_mask, &
             excision_descriptors)
        call tvdreconstruct(nz, ny, nz, xoffset, yoffset, zoffset, &
             Bnz, Bnzplus, Bnzminus, trivial_rp, space_mask, &
             excision_descriptors)
      end if
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if (trivial_rp(i,j,k)) then
              if (flux_direction == 1) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsx, trivialx)
              else if (flux_direction == 2) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsy, trivialy)
              else if (flux_direction == 3) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsz, trivialz)
              end if
            else
              if (flux_direction == 1) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsx, not_trivialx)
              else if (flux_direction == 2) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsy, not_trivialy)
              else if (flux_direction == 3) then
                SpaceMask_SetStateBitsF90(space_mask, i, j, k, \
                                          type_bitsz, not_trivialz)
              end if
            end if
          end do
        end do
      end do
    else
      call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
    end if

  deallocate(trivial_rp)

  where ( (rhoplus < whisky_rho_min).or.(rhominus < whisky_rho_min).or.&
       (epsplus < 0.d0).or.(epsminus < 0.d0) )
      rhoplus = rho
      rhominus = rho
      velxplus = velx
      velxminus = velx
      velyplus = vely
      velyminus = vely
      velzplus = velz
      velzminus = velz
      epsplus = eps
      epsminus = eps
   end where

  if (evolve_tracer .ne. 0) then

    if (use_min_tracer .ne. 0) then
      local_min_tracer = min_tracer
    else
      local_min_tracer = 0d0
    end if
    
    where( (tracerplus  .le. local_min_tracer).or.&
           (tracerminus .le. local_min_tracer) )
      tracerplus = tracer
      tracerminus = tracer
    end where

!!$    write(*,*) 'p2c', local_min_tracer, minval(tracerplus), minval(tracerminus), minval(tracer)

    call CCTK_WARN(0, "evolve_tracer currently broken")
    call Prim2ConservativeTracer(CCTK_PASS_FTOF)

!!$    write(*,*) 'done p2c'

  endif

  return

end subroutine TVD_Reconstruct

