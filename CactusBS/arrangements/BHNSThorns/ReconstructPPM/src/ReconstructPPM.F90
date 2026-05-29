 /*@@
   @file      ReconstructPPM.F90
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
   @routine    Reconstruction
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

subroutine ReconstructPPM(CCTK_ARGUMENTS)
 
  USE PPM_Simple1d_Routines
  USE PPM_Extra_Routines
 
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: nx, ny, nz, i, j, k

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
  !$OMP PARALLEL

  !$OMP WORKSHARE
  lbetax = betax
  lbetay = betay
  lbetaz = betaz
  !$OMP END WORKSHARE NOWAIT

!!$ Initialize variables that store reconstructed quantities

  !$OMP WORKSHARE
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
  !$OMP END WORKSHARE NOWAIT
  if ( whisky_mhd_handle.gt.1 ) then
     !$OMP WORKSHARE
     Bxplus = 0.0d0
     Bxminus = 0.0d0
     Byplus = 0.0d0
     Byminus = 0.0d0
     Bzplus = 0.0d0
     Bzminus = 0.0d0
     !$OMP END WORKSHARE NOWAIT
     if ( clean_divergence .ne. 0 ) then
        !$OMP WORKSHARE
        divclean_psiplus = 0.0d0
        divclean_psiminus = 0.0d0
        !$OMP END WORKSHARE NOWAIT
     end if
  end if

  if (evolve_tracer .ne. 0) then
     !$OMP WORKSHARE
     tracerplus = 0.0d0
     tracerminus = 0.0d0
     !$OMP END WORKSHARE NOWAIT
  endif
  !$OMP END PARALLEL

    if (flux_direction == 1) then
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, nz - whisky_stencil + 1
        do j = whisky_stencil, ny - whisky_stencil + 1
          if (whisky_mhd_handle.gt.1) then

             if (clean_divergence.ne.0) then
                call SimplePPM_1dMHD_wDC(whisky_eos_handle,0,nx,&
                   CCTK_DELTA_SPACE(1),rho(:,j,k),velx(:,j,k),&
                   vely(:,j,k),velz(:,j,k),eps(:,j,k),press(:,j,k),&
                   Bvec(:,j,k,1),Bvec(:,j,k,2),Bvec(:,j,k,3),divclean_psi(:,j,k),&
                   rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                   velzminus(:,j,k),epsminus(:,j,k),Bxminus(:,j,k),&
                   Byminus(:,j,k),Bzminus(:,j,k),divclean_psiminus(:,j,k),&
                   rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),&
                   velzplus(:,j,k),epsplus(:,j,k),Bxplus(:,j,k),&
                   Byplus(:,j,k),Bzplus(:,j,k),divclean_psiplus(:,j,k),&
                   trivial_rp(:,j,k), space_mask(:,j,k), &
                   excision_descriptors,&
                   gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                   gzz(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                   w_lorentz(:,j,k), &
                   1, j, k, nx, ny, nz, mppm_eigenvalue_x_left, &
                                        mppm_eigenvalue_x_right, &
                                        mppm_xwind)
             else
                call SimplePPM_1dMHD(whisky_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                  rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),eps(:,j,k),&
                  press(:,j,k),&
                  Bvec(:,j,k,1),Bvec(:,j,k,2),Bvec(:,j,k,3),&
                  rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                  velzminus(:,j,k),epsminus(:,j,k),&
                  Bxminus(:,j,k),Byminus(:,j,k),Bzminus(:,j,k),&
                  rhoplus(:,j,k),velxplus(:,j,k),velyplus(:,j,k),&
                  velzplus(:,j,k),epsplus(:,j,k),&
                  Bxplus(:,j,k),Byplus(:,j,k),Bzplus(:,j,k),&
                  trivial_rp(:,j,k), space_mask(:,j,k), &
                  excision_descriptors,&
                  gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                  gzz(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                  w_lorentz(:,j,k), &
                  1, j, k, nx, ny, nz, mppm_eigenvalue_x_left, &
                                       mppm_eigenvalue_x_right, &
                                       mppm_xwind)
             end if
          else 
             call SimplePPM_1d(whisky_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                  rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),eps(:,j,k),&
                  press(:,j,k),&
                  rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                  velzminus(:,j,k),epsminus(:,j,k),rhoplus(:,j,k),&
                  velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),epsplus(:,j,k),&
                  trivial_rp(:,j,k), space_mask(:,j,k), &
                  excision_descriptors,&
                  gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                  gzz(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                  w_lorentz(:,j,k), &
                  1, j, k, nx, ny, nz, mppm_eigenvalue_x_left, &
                                       mppm_eigenvalue_x_right, &
                                       mppm_xwind)
          end if
          do i = 1, nx
            if (trivial_rp(i,j,k)) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k,\
                                        type_bitsx, trivialx)
            else
              SpaceMask_SetStateBitsF90(space_mask, i, j, k,\
                                        type_bitsx, not_trivialx)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (j,k)
         do k = whisky_stencil, nz - whisky_stencil + 1
            do j = whisky_stencil, ny - whisky_stencil + 1

               call SimplePPM_tracer_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                    velx(:,j,k),vely(:,j,k),velz(:,j,k), &
                    tracer(:,j,k,:),tracerminus(:,j,k,:),tracerplus(:,j,k,:), & 
                    press(:,j,k))
            end do
         end do
        !$OMP END PARALLEL DO
      end if


    else if (flux_direction == 2) then
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, nz - whisky_stencil + 1
        do j = whisky_stencil, nx - whisky_stencil + 1
          if (whisky_mhd_handle.gt.1) then
             if (clean_divergence.ne.0) then
                call SimplePPM_1dMHD_wDC(whisky_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                  rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),eps(j,:,k),&
                  press(j,:,k),Bvec(j,:,k,2),Bvec(j,:,k,3),Bvec(j,:,k,1),divclean_psi(j,:,k),&
                  rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                  velxminus(j,:,k),epsminus(j,:,k),&
                  Byminus(j,:,k),Bzminus(j,:,k),Bxminus(j,:,k),divclean_psiminus(j,:,k),&
                  rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),&
                  velxplus(j,:,k),epsplus(j,:,k),&
                  Byplus(j,:,k),Bzplus(j,:,k),Bxplus(j,:,k),divclean_psiplus(j,:,k),&
                  trivial_rp(j,:,k), space_mask(j,:,k), &
                  excision_descriptors,&
                  gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                  gxx(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                  w_lorentz(j,:,k), &
                  2, j, k, nx, ny, nz, mppm_eigenvalue_y_left, &
                                       mppm_eigenvalue_y_right, &
                                       mppm_xwind)
             else
               call SimplePPM_1dMHD(whisky_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                  rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),eps(j,:,k),&
                  press(j,:,k),&
                  Bvec(j,:,k,2),Bvec(j,:,k,3),Bvec(j,:,k,1),&
                  rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                  velxminus(j,:,k),epsminus(j,:,k),&
                  Byminus(j,:,k),Bzminus(j,:,k),Bxminus(j,:,k),&
                  rhoplus(j,:,k),velyplus(j,:,k),velzplus(j,:,k),&
                  velxplus(j,:,k),epsplus(j,:,k),&
                  Byplus(j,:,k),Bzplus(j,:,k),Bxplus(j,:,k),&
                  trivial_rp(j,:,k), space_mask(j,:,k), &
                  excision_descriptors,&
                  gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                  gxx(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                  w_lorentz(j,:,k), &
                  2, j, k, nx, ny, nz, mppm_eigenvalue_y_left, &
                                       mppm_eigenvalue_y_right, &
                                       mppm_xwind)
             end if
          else
             call SimplePPM_1d(whisky_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                  rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),eps(j,:,k),&
                  press(j,:,k),&
                  rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                  velxminus(j,:,k),epsminus(j,:,k),rhoplus(j,:,k),&
                  velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),epsplus(j,:,k),&
                  trivial_rp(j,:,k), space_mask(j,:,k), &
                  excision_descriptors,&
                  gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                  gxx(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                  w_lorentz(j,:,k), &
                  2, j, k, nx, ny, nz, mppm_eigenvalue_y_left, &
                                       mppm_eigenvalue_y_right, &
                                       mppm_xwind)
          end if
          do i = 1, ny
            if (trivial_rp(j,i,k)) then
              SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, trivialy)
            else
              SpaceMask_SetStateBitsF90(space_mask, j, i, k, type_bitsy, not_trivialy)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO


      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (j,k)
         do k = whisky_stencil, nz - whisky_stencil + 1
            do j = whisky_stencil, nx - whisky_stencil + 1
               call SimplePPM_tracer_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                    vely(j,:,k),velz(j,:,k),velx(j,:,k), &
                    tracer(j,:,k,:),tracerminus(j,:,k,:),tracerplus(j,:,k,:), & 
                    press(j,:,k))
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

    else if (flux_direction == 3) then
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, ny - whisky_stencil + 1
        do j = whisky_stencil, nx - whisky_stencil + 1
          if (whisky_mhd_handle.gt.1) then
             if (clean_divergence.ne.0) then
               call SimplePPM_1dMHD_wDC(whisky_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                  rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),eps(j,k,:),&
                  press(j,k,:),Bvec(j,k,:,3),Bvec(j,k,:,1),Bvec(j,k,:,2),divclean_psi(j,k,:),&
                  rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
                  velyminus(j,k,:),epsminus(j,k,:),&
                  Bzminus(j,k,:),Bxminus(j,k,:),Byminus(j,k,:),divclean_psiminus(j,k,:),&
                  rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),&
                  velyplus(j,k,:),epsplus(j,k,:),&
                  Bzplus(j,k,:),Bxplus(j,k,:),Byplus(j,k,:),divclean_psiplus(j,k,:),&
                  trivial_rp(j,k,:), space_mask(j,k,:), &
                  excision_descriptors,&
                  gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                  gyy(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
                  w_lorentz(j,k,:), &
                  3, j, k, nx, ny, nz, mppm_eigenvalue_z_left, &
                                       mppm_eigenvalue_z_right, &
                                       mppm_xwind)
             else
                call SimplePPM_1dMHD(whisky_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                  rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),eps(j,k,:),&
                  press(j,k,:),&
                  Bvec(j,k,:,3),Bvec(j,k,:,1),Bvec(j,k,:,2),&
                  rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
                  velyminus(j,k,:),epsminus(j,k,:),&
                  Bzminus(j,k,:),Bxminus(j,k,:),Byminus(j,k,:),&
                  rhoplus(j,k,:),velzplus(j,k,:),velxplus(j,k,:),&
                  velyplus(j,k,:),epsplus(j,k,:),&
                  Bzplus(j,k,:),Bxplus(j,k,:),Byplus(j,k,:),&
                  trivial_rp(j,k,:), space_mask(j,k,:), &
                  excision_descriptors,&
                  gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                  gyy(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
                  w_lorentz(j,k,:), &
                  3, j, k, nx, ny, nz, mppm_eigenvalue_z_left, &
                                       mppm_eigenvalue_z_right, &
                                       mppm_xwind)
             end if
          else
             call SimplePPM_1d(whisky_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
                  rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),eps(j,k,:),&
                  press(j,k,:),&
                  rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
                  velyminus(j,k,:),epsminus(j,k,:),rhoplus(j,k,:),&
                  velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),epsplus(j,k,:),&
                  trivial_rp(j,k,:), space_mask(j,k,:), &
                  excision_descriptors,&
                  gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
                  gyy(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
                  w_lorentz(j,k,:), &
                  3, j, k, nx, ny, nz, mppm_eigenvalue_z_left, &
                                       mppm_eigenvalue_z_right, &
                                       mppm_xwind)
          end if
          do i = 1, nz
            if (trivial_rp(j,k,i)) then
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, trivialz)
            else
              SpaceMask_SetStateBitsF90(space_mask, j, k, i, type_bitsz, not_trivialz)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      if(evolve_tracer.ne.0) then
         !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (j,k)
         do k = whisky_stencil, ny - whisky_stencil + 1
            do j = whisky_stencil, nx - whisky_stencil + 1

               call SimplePPM_tracer_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                    velz(j,k,:),velx(j,k,:),vely(j,k,:), &
                    tracer(j,k,:,:),tracerminus(j,k,:,:),tracerplus(j,k,:,:), & 
                    press(j,k,:))
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

    else
      call CCTK_WARN(0, "Flux direction not x,y,z")
    end if

  deallocate(trivial_rp)
  deallocate(lbetax, lbetay, lbetaz)

  !$OMP PARALLEL WORKSHARE
!!  where ( (rhoplus < whisky_rho_min).or.(rhominus < whisky_rho_min).or.&
!!       (epsplus < 0.d0).or.(epsminus < 0.d0) )
!!      rhoplus = rho
!!      rhominus = rho
!!      velxplus = velx
!!      velxminus = velx
!!      velyplus = vely
!!      velyminus = vely
!!      velzplus = velz
!!      velzminus = velz
!!      epsplus = eps
!!      epsminus = eps
!!  end where
  !$OMP END PARALLEL WORKSHARE

  if (evolve_tracer .ne. 0) then

    if (use_min_tracer .ne. 0) then
      local_min_tracer = min_tracer
    else
      local_min_tracer = 0d0
    end if
    
    !$OMP PARALLEL WORKSHARE
    where( (tracerplus  .le. local_min_tracer).or.&
           (tracerminus .le. local_min_tracer) )
      tracerplus = tracer
      tracerminus = tracer
    end where
    !$OMP END PARALLEL WORKSHARE

!!$    write(*,*) 'p2c', local_min_tracer, minval(tracerplus), minval(tracerminus), minval(tracer)

    call Prim2ConservativeTracer(CCTK_PASS_FTOF)

!!$    write(*,*) 'done p2c'

  endif

  return

end subroutine ReconstructPPM

