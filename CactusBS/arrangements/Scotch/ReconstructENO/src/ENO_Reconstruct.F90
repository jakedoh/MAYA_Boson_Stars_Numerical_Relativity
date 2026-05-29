 /*@@
   @file      ENO_Reconstruct.F90
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
   @routine    ENO_Reconstruct
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

subroutine ENO_Reconstruct(CCTK_ARGUMENTS)
  
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
  CCTK_REAL, dimension(:,:,:),allocatable :: lbetax, lbetay, lbetaz

  CCTK_INT :: ierr

  CCTK_REAL :: local_min_tracer

  allocate(trivial_rp(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Allocation problems with trivial_rp")
  end if
  
  allocate(lbetax(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),&
       lbetay(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),&
       lbetaz(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Allocation problems with lbeta")
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

  lbetax = betax
  lbetay = betay
  lbetaz = betaz

!!$ Initialize variables that store reconstructed quantities

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

  if (evolve_tracer .ne. 0) then
     tracerplus = 0.0d0
     tracerminus = 0.0d0
  endif

     if (evolve_tracer .ne. 0) then
        do itracer=1,number_of_tracers
           call tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
                tracer(:,:,:,itracer), tracerplus(:,:,:,itracer), tracerminus(:,:,:,itracer), trivial_rp, &
                space_mask, excision_descriptors)        
        enddo
     end if

    if (flux_direction == 1) then
      do k = whisky_stencil, cctk_lsh(3) - whisky_stencil + 1
        do j = whisky_stencil, cctk_lsh(2) - whisky_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 rho(:,j,k),rhominus(:,j,k),rhoplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 velx(:,j,k),velxminus(:,j,k),velxplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 vely(:,j,k),velyminus(:,j,k),velyplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 velz(:,j,k),velzminus(:,j,k),velzplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 eps(:,j,k),epsminus(:,j,k),epsplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 dens(:,j,k),densminus(:,j,k),densplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 sx(:,j,k),sxminus(:,j,k),sxplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 sy(:,j,k),syminus(:,j,k),syplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 sz(:,j,k),szminus(:,j,k),szplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(1),&
                 tau(:,j,k),tauminus(:,j,k),tauplus(:,j,k),&
                 trivial_rp(:,j,k), space_mask(:,j,k), &
                 excision_descriptors)
          else
            call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
          end if
          do i = 1, cctk_lsh(1)
            if (trivial_rp(i,j,k)) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, trivialx)
            else
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, not_trivialx)
            end if
          end do
        end do
      end do
    else if (flux_direction == 2) then
      do k = whisky_stencil, cctk_lsh(3) - whisky_stencil + 1
        do j = whisky_stencil, cctk_lsh(1) - whisky_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 rho(j,:,k),rhominus(j,:,k),rhoplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 velx(j,:,k),velxminus(j,:,k),velxplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 vely(j,:,k),velyminus(j,:,k),velyplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 velz(j,:,k),velzminus(j,:,k),velzplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 eps(j,:,k),epsminus(j,:,k),epsplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 dens(j,:,k),densminus(j,:,k),densplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 sx(j,:,k),sxminus(j,:,k),sxplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 sy(j,:,k),syminus(j,:,k),syplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 sz(j,:,k),szminus(j,:,k),szplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(2),&
                 tau(j,:,k),tauminus(j,:,k),tauplus(j,:,k),&
                 trivial_rp(j,:,k), space_mask(j,:,k), &
                 excision_descriptors)
          else
            call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
          end if
          do i = 1, cctk_lsh(2)
            if (trivial_rp(j,i,k)) then
              SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, trivialy)
            else
              SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, not_trivialy)
            end if
          end do
        end do
      end do
    else if (flux_direction == 3) then
      do k = whisky_stencil, cctk_lsh(2) - whisky_stencil + 1
        do j = whisky_stencil, cctk_lsh(1) - whisky_stencil + 1
          if (CCTK_EQUALS(recon_vars,"primitive")) then
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 rho(j,k,:),rhominus(j,k,:),rhoplus(j,k,:),&
                 trivial_rp(j,k,:), space_mask(j,k,:), &
                 excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 velx(j,k,:),velxminus(j,k,:),velxplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 vely(j,k,:),velyminus(j,k,:),velyplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 velz(j,k,:),velzminus(j,k,:),velzplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 eps(j,k,:),epsminus(j,k,:),epsplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
          else if (CCTK_EQUALS(recon_vars,"conservative")) then
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 dens(j,k,:),densminus(j,k,:),densplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 sx(j,k,:),sxminus(j,k,:),sxplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 sy(j,k,:),syminus(j,k,:),syplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 sz(j,k,:),szminus(j,k,:),szplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
            call ENOReconstruct1d(eno_order,cctk_lsh(3),&
                 tau(j,k,:),tauminus(j,k,:),tauplus(j,k,:),&
               trivial_rp(j,k,:), space_mask(j,k,:), &
               excision_descriptors)
          else
            call CCTK_WARN(0, "Variable type to reconstruct not recognized.")
          end if
          do i = 1, cctk_lsh(3)
            if (trivial_rp(j,k,i)) then
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, trivialz)
            else
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, not_trivialz)
            end if
          end do
        end do
      end do
    else
      call CCTK_WARN(0, "Flux direction not x,y,z")
    end if

  deallocate(trivial_rp)
  deallocate(lbetax, lbetay, lbetaz)

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

    call Prim2ConservativeTracer(CCTK_PASS_FTOF)

!!$    write(*,*) 'done p2c'

  endif

  return

end subroutine ENO_Reconstruct

