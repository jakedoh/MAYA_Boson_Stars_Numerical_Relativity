 /*@@
   @file      PPM_fast.F90
   @date      Mon Aug 24 11:32:38 EDT 2009
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
   @routine    FastReconstructPPM
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

subroutine FastReconstructPPM(CCTK_ARGUMENTS)
  
  USE PPM_SimplePPM_1d_XYZ 
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
  CCTK_INT, dimension(3) :: excision_descriptors, bvec_dirs
  CCTK_REAL, dimension(:,:,:),allocatable :: lbetax, lbetay, lbetaz

  CCTK_INT :: ierr

  CCTK_REAL :: local_min_tracer

  if (test_ppm_unroll.ne.0) then
    call ReconstructPPM(CCTK_PASS_FTOF)
    !$OMP PARALLEL WORKSHARE
    rhoplusdiff(:,:,:,flux_direction) = rhoplus
    epsplusdiff(:,:,:,flux_direction) =     epsplus
    velxplusdiff(:,:,:,flux_direction) =     velxplus
    velyplusdiff(:,:,:,flux_direction) =     velyplus
    velzplusdiff(:,:,:,flux_direction) =     velzplus
    rhominusdiff(:,:,:,flux_direction) =     rhominus
    epsminusdiff(:,:,:,flux_direction) =     epsminus
    velxminusdiff(:,:,:,flux_direction) =     velxminus
    velyminusdiff(:,:,:,flux_direction) =     velyminus
    velzminusdiff(:,:,:,flux_direction) =     velzminus
    !$OMP END PARALLEL WORKSHARE
    if ( whisky_mhd_handle.gt.1 ) then
      !$OMP PARALLEL WORKSHARE
      Bxplusdiff(:,:,:,flux_direction)  = Bxplus
      Byplusdiff(:,:,:,flux_direction)  = Byplus
      Bzplusdiff(:,:,:,flux_direction)  = Bzplus
      Bxminusdiff(:,:,:,flux_direction) = Bxminus
      Byminusdiff(:,:,:,flux_direction) = Byminus
      Bzminusdiff(:,:,:,flux_direction) = Bzminus
      !$OMP END PARALLEL WORKSHARE
      if ( clean_divergence.ne.0 ) then
        !$OMP PARALLEL WORKSHARE
        dc_psiplusdiff(:,:,:,flux_direction) = divclean_psiplus
        dc_psiminusdiff(:,:,:,flux_direction) = divclean_psiminus
        !$OMP END PARALLEL WORKSHARE
      end if
    end if
  end if

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

  !$OMP PARALLEL WORKSHARE
  lbetax = betax
  lbetay = betay
  lbetaz = betaz
  !$OMP END PARALLEL WORKSHARE

!!$ Initialize variables that store reconstructed quantities

  !$OMP PARALLEL
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
     bxplus = 0.0d0
     bxminus = 0.0d0
     byplus = 0.0d0
     byminus = 0.0d0
     bzplus = 0.0d0
     bzminus = 0.0d0
     !$OMP END WORKSHARE NOWAIT
     if ( clean_divergence.ne.0 ) then
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
      bvec_dirs = (/ 1,2,3 /)
      if ( whisky_mhd_handle.eq.0 ) then
         call SimplePPM_1d_x(whisky_eos_handle,0,&
             rho,velx,vely,velz,eps,press,&
             rhominus,velxminus,velyminus,velzminus,epsminus,&
             rhoplus,velxplus,velyplus,velzplus,epsplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gxx,gxy,gxz,gyy,gyz,gzz,lbetax,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_x_left, &
                         mppm_eigenvalue_x_right, &
                         mppm_xwind, whisky_stencil)
      else if ( clean_divergence.ne.0 ) then
         call SimplePPM_1d_x_dMHD(whisky_eos_handle,0,&
             rho,velx,vely,velz,eps,press,bvec,bvec_dirs,divclean_psi,&
             rhominus,velxminus,velyminus,velzminus,epsminus,&
             bxminus,byminus,bzminus,divclean_psiminus,&
             rhoplus,velxplus,velyplus,velzplus,epsplus,&
             bxplus,byplus,bzplus,divclean_psiplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gxx,gxy,gxz,gyy,gyz,gzz,lbetax,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_x_left, &
                         mppm_eigenvalue_x_right, &
                         mppm_xwind, whisky_stencil)
      else 
         call SimplePPM_1d_x_MHD(whisky_eos_handle,0,&
             rho,velx,vely,velz,eps,press,bvec,bvec_dirs,&
             rhominus,velxminus,velyminus,velzminus,epsminus,&
             bxminus,byminus,bzminus,&
             rhoplus,velxplus,velyplus,velzplus,epsplus,&
             bxplus,byplus,bzplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gxx,gxy,gxz,gyy,gyz,gzz,lbetax,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_x_left, &
                         mppm_eigenvalue_x_right, &
                         mppm_xwind, whisky_stencil)
      end if
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, nz - whisky_stencil + 1
        do j = whisky_stencil, ny - whisky_stencil + 1
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
         ! Intel 11.1 causes SEGFAULTs when declaring array whose size is unknown
         ! at compile time as FIRSTPRIVATE. We use SHARED insteadd.
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
      bvec_dirs = (/ 2,3,1 /)
      if ( whisky_mhd_handle.eq.0 ) then
         call SimplePPM_1d_y(whisky_eos_handle,0,&
             rho,vely,velz,velx,eps,press,&
             rhominus,velyminus,velzminus,velxminus,epsminus,&
             rhoplus,velyplus,velzplus,velxplus,epsplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gyy,gyz,gxy,gzz,gxz,gxx,lbetay,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_y_left, &
                         mppm_eigenvalue_y_right, &
                         mppm_xwind, whisky_stencil)
      else if ( clean_divergence.ne.0 ) then
         call SimplePPM_1d_y_dMHD(whisky_eos_handle,0,&
             rho,vely,velz,velx,eps,press,bvec,bvec_dirs,divclean_psi,&
             rhominus,velyminus,velzminus,velxminus,epsminus,&
             byminus,bzminus,bxminus,divclean_psiminus,&
             rhoplus,velyplus,velzplus,velxplus,epsplus,&
             byplus,bzplus,bxplus,divclean_psiplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gyy,gyz,gxy,gzz,gxz,gxx,lbetay,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_y_left, &
                         mppm_eigenvalue_y_right, &
                         mppm_xwind, whisky_stencil)
      else 
         call SimplePPM_1d_y_MHD(whisky_eos_handle,0,&
             rho,vely,velz,velx,eps,press,bvec,bvec_dirs,&
             rhominus,velyminus,velzminus,velxminus,epsminus,&
             byminus,bzminus,bxminus,&
             rhoplus,velyplus,velzplus,velxplus,epsplus,&
             byplus,bzplus,bxplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gyy,gyz,gxy,gzz,gxz,gxx,lbetay,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_y_left, &
                         mppm_eigenvalue_y_right, &
                         mppm_xwind, whisky_stencil)
      end if
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, nz - whisky_stencil + 1
        do j = whisky_stencil, nx - whisky_stencil + 1
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
      bvec_dirs = (/ 3,1,2 /)
      if ( whisky_mhd_handle.eq.0 ) then
         call SimplePPM_1d_z(whisky_eos_handle,0,&
             rho,velz,velx,vely,eps,press,&
             rhominus,velzminus,velxminus,velyminus,epsminus,&
             rhoplus,velzplus,velxplus,velyplus,epsplus,&
             trivial_rp,space_mask,excision_descriptors,&
             gzz,gxz,gyz,gxx,gxy,gyy,lbetaz,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_z_left, &
                         mppm_eigenvalue_z_right, &
                         mppm_xwind, whisky_stencil)
      else if ( clean_divergence.ne.0 ) then
         call SimplePPM_1d_z_dMHD(whisky_eos_handle,0,&
             rho,velz,velx,vely,eps,press,bvec,bvec_dirs,divclean_psi,&
             rhominus,velzminus,velxminus,velyminus,epsminus,&
             bzminus,bxminus,byminus,divclean_psiminus,&
             rhoplus,velzplus,velxplus,velyplus,epsplus,&
             bzplus,bxplus,byplus,divclean_psiplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gzz,gxz,gyz,gxx,gxy,gyy,lbetaz,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_z_left, &
                         mppm_eigenvalue_z_right, &
                         mppm_xwind, whisky_stencil)
      else 
         call SimplePPM_1d_z_MHD(whisky_eos_handle,0,&
             rho,velz,velx,vely,eps,press,bvec,bvec_dirs,&
             rhominus,velzminus,velxminus,velyminus,epsminus,&
             bzminus,bxminus,byminus,&
             rhoplus,velzplus,velxplus,velyplus,epsplus,&
             bzplus,bxplus,byplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gzz,gxz,gyz,gxx,gxy,gyy,lbetaz,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_z_left, &
                         mppm_eigenvalue_z_right, &
                         mppm_xwind, whisky_stencil)
      end if
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, ny - whisky_stencil + 1
        do j = whisky_stencil, nx - whisky_stencil + 1
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

!!  !$OMP PARALLEL WORKSHARE
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
!!  !$OMP END PARALLEL WORKSHARE

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

  if (test_ppm_unroll.ne.0) then
    !$OMP PARALLEL WORKSHARE
    rhoplusdiff(:,:,:,flux_direction)   = rhoplusdiff(:,:,:,flux_direction)   - rhoplus
    epsplusdiff(:,:,:,flux_direction)   = epsplusdiff(:,:,:,flux_direction)   - epsplus
    velxplusdiff(:,:,:,flux_direction)  = velxplusdiff(:,:,:,flux_direction)  - velxplus
    velyplusdiff(:,:,:,flux_direction)  = velyplusdiff(:,:,:,flux_direction)  - velyplus
    velzplusdiff(:,:,:,flux_direction)  = velzplusdiff(:,:,:,flux_direction)  - velzplus
    rhominusdiff(:,:,:,flux_direction)  = rhominusdiff(:,:,:,flux_direction)  - rhominus
    epsminusdiff(:,:,:,flux_direction)  = epsminusdiff(:,:,:,flux_direction)  - epsminus
    velxminusdiff(:,:,:,flux_direction) = velxminusdiff(:,:,:,flux_direction) - velxminus
    velyminusdiff(:,:,:,flux_direction) = velyminusdiff(:,:,:,flux_direction) - velyminus
    velzminusdiff(:,:,:,flux_direction) = velzminusdiff(:,:,:,flux_direction) - velzminus
    !$OMP END PARALLEL WORKSHARE
    if ( whisky_mhd_handle.gt.1 ) then
      !$OMP PARALLEL WORKSHARE
      Bxplusdiff(:,:,:,flux_direction)  = Bxplusdiff(:,:,:,flux_direction)    - Bxplus
      Byplusdiff(:,:,:,flux_direction)  = Byplusdiff(:,:,:,flux_direction)    - Byplus
      Bzplusdiff(:,:,:,flux_direction)  = Bzplusdiff(:,:,:,flux_direction)    - Bzplus
      Bxminusdiff(:,:,:,flux_direction) = Bxminusdiff(:,:,:,flux_direction)   - Bxminus
      Byminusdiff(:,:,:,flux_direction) = Byminusdiff(:,:,:,flux_direction)   - Byminus
      Bzminusdiff(:,:,:,flux_direction) = Bzminusdiff(:,:,:,flux_direction)   - Bzminus
      !$OMP END PARALLEL WORKSHARE
      if ( clean_divergence.ne.0 ) then
        !$OMP PARALLEL WORKSHARE
        dc_psiplusdiff(:,:,:,flux_direction) = dc_psiplusdiff(:,:,:,flux_direction)   - divclean_psiplus
        dc_psiminusdiff(:,:,:,flux_direction) = dc_psiminusdiff(:,:,:,flux_direction) - divclean_psiminus
        !$OMP END PARALLEL WORKSHARE
      end if
    end if
  end if

  return

end subroutine FastReconstructPPM

 /*@@
   @routine    FastPPM_Polytype
   @date       Tue Mar 19 11:36:55 2002
   @author     Ian Hawke  
   @desc 
   If using a polytropic type EOS, we do not have to reconstruct all the 
   variables.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine FastPPM_Polytype(CCTK_ARGUMENTS)
 
  USE PPM_SimplePPM_1d_XYZ
  USE PPM_Extra_Routines
 
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: nx, ny, nz, i, j, k

  logical, dimension(:,:,:), allocatable :: trivial_rp
!!$  logical, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: trivial_rp

  CCTK_INT :: type_bitsx, trivialx, not_trivialx, &
       &type_bitsy, trivialy, not_trivialy, &
       &type_bitsz, trivialz, not_trivialz

  CCTK_INT, dimension(3) :: excision_descriptors, bvec_dirs
  CCTK_REAL, dimension(:,:,:),allocatable :: lbetax, lbetay, lbetaz

  CCTK_INT :: ierr

  if (test_ppm_unroll.ne.0) then
    call PPM_Polytype(CCTK_PASS_FTOF)
    !$OMP PARALLEL WORKSHARE
    rhoplusdiff(:,:,:,flux_direction) = rhoplus
    epsplusdiff(:,:,:,flux_direction) =     epsplus
    velxplusdiff(:,:,:,flux_direction) =     velxplus
    velyplusdiff(:,:,:,flux_direction) =     velyplus
    velzplusdiff(:,:,:,flux_direction) =     velzplus
    rhominusdiff(:,:,:,flux_direction) =     rhominus
    epsminusdiff(:,:,:,flux_direction) =     epsminus
    velxminusdiff(:,:,:,flux_direction) =     velxminus
    velyminusdiff(:,:,:,flux_direction) =     velyminus
    velzminusdiff(:,:,:,flux_direction) =     velzminus
    !$OMP END PARALLEL WORKSHARE
    if ( whisky_mhd_handle.gt.1 ) then
      !$OMP PARALLEL WORKSHARE
      Bxplusdiff(:,:,:,flux_direction)  = Bxplus
      Byplusdiff(:,:,:,flux_direction)  = Byplus
      Bzplusdiff(:,:,:,flux_direction)  = Bzplus
      Bxminusdiff(:,:,:,flux_direction) = Bxminus
      Byminusdiff(:,:,:,flux_direction) = Byminus
      Bzminusdiff(:,:,:,flux_direction) = Bzminus
      !$OMP END PARALLEL WORKSHARE
      if ( clean_divergence.ne.0 ) then
        !$OMP PARALLEL WORKSHARE
        dc_psiplusdiff(:,:,:,flux_direction) = divclean_psiplus
        dc_psiminusdiff(:,:,:,flux_direction) = divclean_psiminus
        !$OMP END PARALLEL WORKSHARE
      end if
    end if
  end if

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

  !$OMP PARALLEL
  !$OMP WORKSHARE
  trivial_rp = .false.
  !$OMP END WORKSHARE NOWAIT

!!$  Currently only option is reconstruction on primitive variables.
!!$  Should change this.

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
     bxplus = 0.0d0
     bxminus = 0.0d0
     byplus = 0.0d0
     byminus = 0.0d0
     bzplus = 0.0d0
     bzminus = 0.0d0
     !$OMP END WORKSHARE NOWAIT
  end if
  if ( clean_divergence.ne.0 ) then
     !$OMP WORKSHARE
     divclean_psiplus = 0.0d0
     divclean_psiminus = 0.0d0
     !$OMP END WORKSHARE NOWAIT
  end if

  if (evolve_tracer .ne. 0) then
     !$OMP WORKSHARE
     tracerplus = 0.0d0
     tracerminus = 0.0d0
     !$OMP END WORKSHARE NOWAIT
  end if
  !$OMP END PARALLEL

    if (flux_direction == 1) then
      bvec_dirs = (/ 1,2,3 /)
      if ( whisky_mhd_handle.eq.0 ) then
         call SimplePPM_1d_x(whisky_eos_handle,1,&
             rho,velx,vely,velz,eps,press,&
             rhominus,velxminus,velyminus,velzminus,epsminus,&
             rhoplus,velxplus,velyplus,velzplus,epsplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gxx,gxy,gxz,gyy,gyz,gzz,lbetax,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_x_left, &
                         mppm_eigenvalue_x_right, &
                         mppm_xwind, whisky_stencil)
      else if ( clean_divergence.ne.0 ) then
         call SimplePPM_1d_x_dMHD(whisky_eos_handle,1,&
             rho,velx,vely,velz,eps,press,bvec,bvec_dirs,divclean_psi,&
             rhominus,velxminus,velyminus,velzminus,epsminus,&
             bxminus,byminus,bzminus,divclean_psiminus,&
             rhoplus,velxplus,velyplus,velzplus,epsplus,&
             bxplus,byplus,bzplus,divclean_psiplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gxx,gxy,gxz,gyy,gyz,gzz,lbetax,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_x_left, &
                         mppm_eigenvalue_x_right, &
                         mppm_xwind, whisky_stencil)
      else 
         call SimplePPM_1d_x_MHD(whisky_eos_handle,1,&
             rho,velx,vely,velz,eps,press,bvec,bvec_dirs,&
             rhominus,velxminus,velyminus,velzminus,epsminus,&
             bxminus,byminus,bzminus,&
             rhoplus,velxplus,velyplus,velzplus,epsplus,&
             bxplus,byplus,bzplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gxx,gxy,gxz,gyy,gyz,gzz,lbetax,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_x_left, &
                         mppm_eigenvalue_x_right, &
                         mppm_xwind, whisky_stencil)
      end if
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, nz - whisky_stencil + 1
        do j = whisky_stencil, ny - whisky_stencil + 1

          if(evolve_tracer.ne.0) then
             call SimplePPM_tracer_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                  velx(:,j,k),vely(:,j,k),velz(:,j,k), &
                  tracer(:,j,k,:),tracerminus(:,j,k,:),tracerplus(:,j,k,:), &
                  press(:,j,k))
          end if

          do i = 1, nx
            if (trivial_rp(i,j,k)) then
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, trivialx)
            else
              SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bitsx, not_trivialx)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    else if (flux_direction == 2) then
      bvec_dirs = (/ 2,3,1 /)
      if ( whisky_mhd_handle.eq.0 ) then
         call SimplePPM_1d_y(whisky_eos_handle,1,&
             rho,vely,velz,velx,eps,press,&
             rhominus,velyminus,velzminus,velxminus,epsminus,&
             rhoplus,velyplus,velzplus,velxplus,epsplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gyy,gyz,gxy,gzz,gxz,gxx,lbetay,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_y_left, &
                         mppm_eigenvalue_y_right, &
                         mppm_xwind, whisky_stencil)
      else if ( clean_divergence.ne.0 ) then
         call SimplePPM_1d_y_dMHD(whisky_eos_handle,1,&
             rho,vely,velz,velx,eps,press,bvec,bvec_dirs,divclean_psi,&
             rhominus,velyminus,velzminus,velxminus,epsminus,&
             byminus,bzminus,bxminus,divclean_psiminus,&
             rhoplus,velyplus,velzplus,velxplus,epsplus,&
             byplus,bzplus,bxplus,divclean_psiplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gyy,gyz,gxy,gzz,gxz,gxx,lbetay,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_y_left, &
                         mppm_eigenvalue_y_right, &
                         mppm_xwind, whisky_stencil)
      else 
         call SimplePPM_1d_y_MHD(whisky_eos_handle,1,&
             rho,vely,velz,velx,eps,press,bvec,bvec_dirs,&
             rhominus,velyminus,velzminus,velxminus,epsminus,&
             byminus,bzminus,bxminus,&
             rhoplus,velyplus,velzplus,velxplus,epsplus,&
             byplus,bzplus,bxplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gyy,gyz,gxy,gzz,gxz,gxx,lbetay,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_y_left, &
                         mppm_eigenvalue_y_right, &
                         mppm_xwind, whisky_stencil)
      end if
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, nz - whisky_stencil + 1
        do j = whisky_stencil, nx - whisky_stencil + 1

          if(evolve_tracer.ne.0) then
             call SimplePPM_tracer_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                  vely(j,:,k),velz(j,:,k),velx(j,:,k), &
                  tracer(j,:,k,:),tracerminus(j,:,k,:),tracerplus(j,:,k,:), &
                  press(j,:,k))
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
    else if (flux_direction == 3) then
      bvec_dirs = (/ 3,1,2 /)
      if ( whisky_mhd_handle.eq.0 ) then
         call SimplePPM_1d_z(whisky_eos_handle,1,&
             rho,velz,velx,vely,eps,press,&
             rhominus,velzminus,velxminus,velyminus,epsminus,&
             rhoplus,velzplus,velxplus,velyplus,epsplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gzz,gxz,gyz,gxx,gxy,gyy,lbetaz,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_z_left, &
                         mppm_eigenvalue_z_right, &
                         mppm_xwind, whisky_stencil)
      else if ( clean_divergence.ne.0 ) then
         call SimplePPM_1d_z_dMHD(whisky_eos_handle,1,&
             rho,velz,velx,vely,eps,press,bvec,bvec_dirs,divclean_psi,&
             rhominus,velzminus,velxminus,velyminus,epsminus,&
             bzminus,bxminus,byminus,divclean_psiminus,&
             rhoplus,velzplus,velxplus,velyplus,epsplus,&
             bzplus,bxplus,byplus,divclean_psiplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gzz,gxz,gyz,gxx,gxy,gyy,lbetaz,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_z_left, &
                         mppm_eigenvalue_z_right, &
                         mppm_xwind, whisky_stencil)
      else 
         call SimplePPM_1d_z_MHD(whisky_eos_handle,1,&
             rho,velz,velx,vely,eps,press,bvec,bvec_dirs,&
             rhominus,velzminus,velxminus,velyminus,epsminus,&
             bzminus,bxminus,byminus,&
             rhoplus,velzplus,velxplus,velyplus,epsplus,&
             bzplus,bxplus,byplus,&
             trivial_rp,space_mask, excision_descriptors,&
             gzz,gxz,gyz,gxx,gxy,gyy,lbetaz,alp,w_lorentz,&
             nx, ny, nz, mppm_eigenvalue_z_left, &
                         mppm_eigenvalue_z_right, &
                         mppm_xwind, whisky_stencil)
      end if
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE (i,j,k)
      do k = whisky_stencil, ny - whisky_stencil + 1
        do j = whisky_stencil, nx - whisky_stencil + 1
          
          if(evolve_tracer.ne.0) then
             call SimplePPM_tracer_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                  velz(j,k,:),velx(j,k,:),vely(j,k,:), &
                  tracer(j,k,:,:),tracerminus(j,k,:,:),tracerplus(j,k,:,:), &
                  press(j,k,:))
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
    else
      call CCTK_WARN(0, "Flux direction not x,y,z")
    end if

  deallocate(trivial_rp)
  deallocate(lbetax, lbetay, lbetaz)

!!  !$OMP PARALLEL WORKSHARE
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
!!  !$OMP END PARALLEL WORKSHARE

    if (evolve_tracer .ne. 0) then
      call Prim2ConservativeTracer(CCTK_PASS_FTOF)
    end if

  if (test_ppm_unroll.ne.0) then
    !$OMP PARALLEL WORKSHARE
    rhoplusdiff(:,:,:,flux_direction)   = rhoplusdiff(:,:,:,flux_direction)   - rhoplus
    epsplusdiff(:,:,:,flux_direction)   = epsplusdiff(:,:,:,flux_direction)   - epsplus
    velxplusdiff(:,:,:,flux_direction)  = velxplusdiff(:,:,:,flux_direction)  - velxplus
    velyplusdiff(:,:,:,flux_direction)  = velyplusdiff(:,:,:,flux_direction)  - velyplus
    velzplusdiff(:,:,:,flux_direction)  = velzplusdiff(:,:,:,flux_direction)  - velzplus
    rhominusdiff(:,:,:,flux_direction)  = rhominusdiff(:,:,:,flux_direction)  - rhominus
    epsminusdiff(:,:,:,flux_direction)  = epsminusdiff(:,:,:,flux_direction)  - epsminus
    velxminusdiff(:,:,:,flux_direction) = velxminusdiff(:,:,:,flux_direction) - velxminus
    velyminusdiff(:,:,:,flux_direction) = velyminusdiff(:,:,:,flux_direction) - velyminus
    velzminusdiff(:,:,:,flux_direction) = velzminusdiff(:,:,:,flux_direction) - velzminus
    !$OMP END PARALLEL WORKSHARE
    if ( whisky_mhd_handle.gt.1 ) then
      !$OMP PARALLEL WORKSHARE
      Bxplusdiff(:,:,:,flux_direction)  = Bxplusdiff(:,:,:,flux_direction)    - Bxplus
      Byplusdiff(:,:,:,flux_direction)  = Byplusdiff(:,:,:,flux_direction)    - Byplus
      Bzplusdiff(:,:,:,flux_direction)  = Bzplusdiff(:,:,:,flux_direction)    - Bzplus
      Bxminusdiff(:,:,:,flux_direction) = Bxminusdiff(:,:,:,flux_direction)   - Bxminus
      Byminusdiff(:,:,:,flux_direction) = Byminusdiff(:,:,:,flux_direction)   - Byminus
      Bzminusdiff(:,:,:,flux_direction) = Bzminusdiff(:,:,:,flux_direction)   - Bzminus
      !$OMP END PARALLEL WORKSHARE
      if ( clean_divergence.ne.0 ) then
        !$OMP PARALLEL WORKSHARE
        dc_psiplusdiff(:,:,:,flux_direction) = dc_psiplusdiff(:,:,:,flux_direction)   - divclean_psiplus
        dc_psiminusdiff(:,:,:,flux_direction) = dc_psiminusdiff(:,:,:,flux_direction) - divclean_psiminus
        !$OMP END PARALLEL WORKSHARE
      end if
    end if
  end if

  return

end subroutine FastPPM_Polytype
